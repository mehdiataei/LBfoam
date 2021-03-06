/*
 *  LBfoam: An open-source software package for the simulation of foaming
 *  using the Lattice Boltzmann Method
 *  Copyright (C) 2020 Mohammadmehdi Ataei
 *  m.ataei@mail.utoronto.ca
 *  This file is part of LBfoam.
 *
 *  LBfoam is free software: you can redistribute it and/or modify it under
 *  the terms of the GNU Affero General Public License as published by the
 *  Free Software Foundation version 3.
 *
 *  LBfoam is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
 *  more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this Program. If not, see <http://www.gnu.org/licenses/>.
 *
 *  #############################################################################
 *
 *  Author:         Mohammadmehdi Ataei, 2020
 *
 *  #############################################################################
 *
 *  Parts of the LBfoam code that originate from Palabos are distributed
 *  under the terms of the AGPL 3.0 license with the following copyright
 *  notice:
 *
 *  This file is part of the Palabos library.
 *
 *  Copyright (C) 2011-2017 FlowKit Sarl
 *  Route d'Oron 2
 *  1010 Lausanne, Switzerland
 *  E-mail contact: contact@flowkit.com
 *
 *  The most recent release of Palabos can be downloaded at
 *  <http://www.palabos.org/>
 *
 *  The library Palabos is free software: you can redistribute it and/or
 *  modify it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  The library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "lbfoam/algorithms/PLIC2D.h"
#include "lbfoam/algorithms/PLIC3D.h"
#include "lbfoam/algorithms/poisson_disk_sampling.h"
#include "lbfoam/algorithms/rayTracer2D.h"
#include "lbfoam/bubble/bubbleGrowth2D.h"
#include "lbfoam/bubble/bubbleTracking2D.h"
#include "lbfoam/dynamics/advectionDiffusionBubbleGrowth2D.h"
#include "lbfoam/models/bodyForce2D.h"
#include "lbfoam/models/createBubbles2D.h"
#include "lbfoam/models/freeSurfaceAnalysis2D.h"
#include "lbfoam/models/freeSurfaceBoundaryCondition2D.h"
#include "lbfoam/models/freeSurfaceInitializer2D.h"
#include "lbfoam/models/freeSurfaceModel2D.h"
#include "lbfoam/models/freeSurfaceTemplates2D.h"
#include "lbfoam/models/freeSurfaceUtil2D.h"
#include "lbfoam/models/immersedWalls2D.h"