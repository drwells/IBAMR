// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files
#include <IBAMR_config.h>
// #include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// other samrai stuff
#include <HierarchyDataOpsManager.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_partitioner.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/BoxPartitioner.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// This file is the main driver for force spreading tests (i.e.,
// IBFEmethod::spreadForce). At the moment it simply prints out the force
// values.

// Coordinate mapping function.
void
coordinate_mapping_function(libMesh::Point& X, const libMesh::Point& s, void* /*ctx*/)
{
    X(0) = s(0) + 0.6;
    X(1) = s(1) + 0.5;
#if (NDIM == 3)
    X(2) = s(2) + 0.5;
#endif
    return;
} // coordinate_mapping_function

int
main(int argc, char** argv)
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-16");
    PetscOptionsSetValue(nullptr, "-ksp_atol", "1e-16");

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create a simple FE mesh.
        ReplicatedMesh mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const double R = 0.2;
        if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
        }
        mesh.prepare_for_use();
        // metis does a good job partitioning, but the partitioning relies on
        // random numbers: the seed changed in libMesh commit
        // 98cede90ca8837688ee13aac5e299a3765f083da (between 1.3.1 and
        // 1.4.0). Hence, to achieve consistent partitioning, use a simpler partitioning scheme:
        LinearPartitioner partitioner;
        partitioner.partition(mesh);

#if 0
        {
            std::cout.precision(16);
            auto node_begin = mesh.active_nodes_begin();
            auto node_end = mesh.active_nodes_end();
            for (auto n_it = node_begin; n_it != node_end; ++n_it)
                std::cout << static_cast<TypeVector<double> &>(**n_it)(0) << ", "
                          << static_cast<TypeVector<double> &>(**n_it)(1) << '\n';
        }
#endif

        plog << "Number of elements: " << mesh.n_active_elem() << std::endl;

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"), false);
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry, false);
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"),
                false);
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"),
                false);
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           false);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator,
                                              false);

        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer,
                                        false);

        // Configure the IBFE solver.
        ib_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);
        ib_method_ops->initializeFEEquationSystems();
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // only regrid if we are in parallel - do this to avoid weird problems
        // with inconsistent initial partitionings (see
        // IBFEMethod::d_skip_initial_workload_log)
        if (SAMRAI_MPI::getNodes() != 1) time_integrator->regridHierarchy();

        // Now for the actual test. Set up a new variable containing ghost data:
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const Pointer<SAMRAI::hier::Variable<NDIM> > f_var = time_integrator->getBodyForceVariable();
        const Pointer<VariableContext> f_ghost_ctx = var_db->getContext("f_ghost");

        int n_ghosts = 3;
        if (app_initializer->getComponentDatabase("IBFEMethod")->keyExists("min_ghost_cell_width"))
        {
            n_ghosts = app_initializer->getComponentDatabase("IBFEMethod")->getInteger("min_ghost_cell_width");
        }
        const int f_ghost_idx = var_db->registerVariableAndContext(f_var, f_ghost_ctx, n_ghosts);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(f_ghost_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_ghost_idx);
                f_data->fillAll(0.0);
            }
        }

        // synch ghost data
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
        ghost_cell_components[0] = InterpolationTransactionComponent(f_ghost_idx,
                                                                     "CONSERVATIVE_LINEAR_REFINE",
                                                                     true,
                                                                     "CONSERVATIVE_COARSEN",
                                                                     "LINEAR",
                                                                     false,
                                                                     {}, // f_bc_coefs
                                                                     nullptr);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/ 0.0);

        const double dt = time_integrator->getMaximumTimeStepSize();
        time_integrator->preprocessIntegrateHierarchy(
            time_integrator->getIntegratorTime(), time_integrator->getIntegratorTime() + dt, 1 /*???*/);
        auto& fe_data_manager = *ib_method_ops->getFEDataManager();
        auto& equation_systems = *fe_data_manager.getEquationSystems();
        auto& force_system = equation_systems.get_system(IBFEMethod::FORCE_SYSTEM_NAME);
        auto& half_f_vector = dynamic_cast<libMesh::PetscVector<double>&>(*force_system.current_local_solution);
        for (unsigned int i = half_f_vector.first_local_index(); i < half_f_vector.last_local_index(); ++i)
        {
            half_f_vector.set(i, i % 10);
        }
        half_f_vector.close();

        std::ostringstream out;
        if (SAMRAI_MPI::getNodes() != 1)
        {
            // partitioning is only relevant when there are multiple processors
            Pointer<PatchLevel<NDIM> > patch_level =
                patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
            const BoxArray<NDIM> boxes = patch_level->getBoxes();
            plog << "hierarchy boxes:\n";
            for (int i = 0; i < boxes.size(); ++i) plog << boxes[i] << '\n';
            // rank is only relevant when there are multiple processors
            out << "\nrank: " << SAMRAI_MPI::getRank() << '\n';
        }

        ib_method_ops->spreadForce(f_ghost_idx, nullptr, {}, time_integrator->getIntegratorTime() + dt / 2);
        const double cutoff = input_db->getDoubleWithDefault("output_cutoff_value", 0.0);
        {
            const int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                bool printed_value = false;
                std::ostringstream patch_out;
                patch_out << "patch number " << p() << '\n';
                patch_out.precision(16);
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                patch_out << patch->getBox() << '\n';
                Pointer<CartesianPatchGeometry<NDIM> > patch_geo = patch->getPatchGeometry();
                patch_out << "patch x_lo: " << patch_geo->getXLower()[0] << ", "
                          << patch_geo->getXLower()[1] << '\n';
                patch_out << "patch x_up: " << patch_geo->getXUpper()[0] << ", "
                          << patch_geo->getXUpper()[1] << '\n';

                Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_ghost_idx);
                const Box<NDIM> patch_box = patch->getBox();

                // same as SideData::print, but elides zero values. We don't
                // print any information about the patch when no values are
                // above the cutoff.
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    patch_out << "Array side normal = " << axis << std::endl;
                    for (int d = 0; d < f_data->getDepth(); ++d)
                    {
                        patch_out << "Array depth = " << d << std::endl;
                        const ArrayData<NDIM, double>& data = f_data->getArrayData(axis);
                        for (SideIterator<NDIM> i(patch_box, axis); i; i++)
                        {
                            const double value = data(i(), d);
                            if (std::abs(value) > cutoff)
                            {
                                patch_out << "array" << i() << " = " << value << '\n';
                                printed_value = true;
                            }
                        }
                    }
                }
                if (printed_value) out << patch_out.str();
                // f_data->print(patch_box, plog, 16);
            }
        }
        SAMRAI_MPI::barrier();

        // okay, now each processor has some output, but we want to get
        // everything on rank 0 to print to plog.
        const std::string to_log = out.str();
        const int n_nodes = SAMRAI_MPI::getNodes();
        std::vector<unsigned long> string_sizes(n_nodes);

        const unsigned long size = to_log.size();
        int ierr = MPI_Gather(
            &size, 1, MPI_UNSIGNED_LONG, string_sizes.data(), 1, MPI_UNSIGNED_LONG, 0, SAMRAI_MPI::getCommunicator());
        TBOX_ASSERT(ierr == 0);

        // MPI_Gatherv would be more efficient, but this just a test so its
        // not too important
        if (SAMRAI_MPI::getRank() == 0)
        {
            plog << to_log;
            for (int r = 1; r < n_nodes; ++r)
            {
                std::string input;
                input.resize(string_sizes[r]);
                ierr = MPI_Recv(
                    &input[0], string_sizes[r], MPI_CHAR, r, 0, SAMRAI_MPI::getCommunicator(), MPI_STATUS_IGNORE);
                TBOX_ASSERT(ierr == 0);
                plog << input;
            }
        }
        else
            MPI_Send(to_log.data(), size, MPI_CHAR, 0, 0, SAMRAI_MPI::getCommunicator());
    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
} // main
