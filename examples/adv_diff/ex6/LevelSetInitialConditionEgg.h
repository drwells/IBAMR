// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2018 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_LevelSetInitialConditionEgg
#define included_LevelSetInitialConditionEgg

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class LevelSetInitialConditionEgg provides an initial condition for
 * the level set function.
 */
class LevelSetInitialConditionEgg : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LevelSetInitialConditionEgg(const std::string& object_name, const IBTK::VectorNd& origin);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetInitialConditionEgg() = default;

    /*!
     * \brief Indicates whether the concrete LevelSetInitialConditionEgg object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Indicates whether the LevelSetInitialConditionEgg corresponds to inner
     * or outer cylinder.
     */
    void setCylinderType(std::string type);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL)) override;

    //\}

private:
    /*!
     * Deleted default constructor.
     */
    LevelSetInitialConditionEgg() = delete;

    /*!
     * Deleted copy constructor.
     */
    LevelSetInitialConditionEgg(const LevelSetInitialConditionEgg& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    LevelSetInitialConditionEgg& operator=(const LevelSetInitialConditionEgg& that) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Origin of geometry.
     */
    IBTK::VectorNd d_origin;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LevelSetInitialConditionEgg
