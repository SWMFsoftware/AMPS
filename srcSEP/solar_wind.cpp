
#include "sep.h"

namespace SEP {
    namespace SolarWind {
        /**
         * Offset for solar wind velocity divergence in the data buffer
         */
        int DivSolarWindVelocityOffset = -1;

        /**
         * Calculates the solar wind velocity at a given position.
         * 
         * @param v Output array for velocity components
         * @param x Position vector
         * @return The magnitude of the solar wind velocity
         */
        double GetSolarWindVelocity(double *v, double *x) {
            // Calculate the length of the position vector
            double length = Vector3D::Length(x);
            
            // Set each component of velocity
            // v[idim] = 400 km/s * x[idim] / |x|
            for (int idim = 0; idim < 3; idim++) {
                v[idim] = 400E3 * x[idim] / length;
            }
            
            // Return the magnitude of the velocity (which is 400 km/s)
            return 400E3;
        }

	/*
         * Request data buffer for storing div(Solar Wind Velocity) 
         * 
         * @param offset The current offset in the data buffer
         * @return The size of allocated buffer in bytes
         */
        int RequestDataBuffer(int offset) {
            DivSolarWindVelocityOffset = offset;
            int TotalDataLength = 1;
            
            return TotalDataLength * sizeof(double);
        }

        void Init() {
            //request sampling buffer and particle fields
            PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestDataBuffer);
        }

     /**
         * Calculates the divergence of the solar wind velocity at a given position
         * 
         * @param x Position vector
         * @param node Pointer to the AMR tree node
         * @return The divergence of the solar wind velocity
         */
        double GetDivSolarWindVelocity(double *x, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
            // Calculate cell size from node dimensions and number of cells
            double cellSize[3];
            cellSize[0] = (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_;
            cellSize[1] = (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_;
            cellSize[2] = (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_;
            
            // Define step size as half of the local cell size for each dimension
            double h[3];
            for (int i = 0; i < 3; i++) {
                h[i] = 0.5 * cellSize[i];
            }
            
            // Calculate the divergence using central differences
            double div = 0.0;
            double vPlus[3], vMinus[3];
            double xPlus[3], xMinus[3];
            
            for (int dim = 0; dim < 3; dim++) {
                // Copy the original position
                for (int i = 0; i < 3; i++) {
                    xPlus[i] = x[i];
                    xMinus[i] = x[i];
                }
                
                // Shift position forward and backward in current dimension
                xPlus[dim] += h[dim];
                xMinus[dim] -= h[dim];
                
                // Check if points are valid (not outside domain or inside Sun)
                bool validPlus = IsPointValid(xPlus);
                bool validMinus = IsPointValid(xMinus);
                
                if (!validPlus || !validMinus) {
                    // If either point is invalid, skip this dimension
                    continue;
                }
                
                // Get velocity at shifted positions
                GetSolarWindVelocity(vPlus,xPlus);
                GetSolarWindVelocity(vMinus, xMinus);
                
                // Calculate partial derivative: dv_i/dx_i using central difference
                div += (vPlus[dim] - vMinus[dim]) / (2.0 * h[dim]);
            }
            
            return div;
        }
        
        /**
         * Check if a point is valid (not outside domain or inside the Sun)
         * 
         * @param x Position vector to check
         * @return true if valid, false otherwise
         */
        bool IsPointValid(double *x) {
            if (Vector3D::DotProduct(x,x) < _SUN__RADIUS_*_SUN__RADIUS_) {
                return false;
            }
            
            // Check if outside domain
            // Use the global domain boundaries defined in PIC
            const double *xmin = PIC::Mesh::mesh->xGlobalMin;
            const double *xmax = PIC::Mesh::mesh->xGlobalMax;
            
            if (x[0] < xmin[0] || x[0] > xmax[0] ||
                x[1] < xmin[1] || x[1] > xmax[1] ||
                x[2] < xmin[2] || x[2] > xmax[2]) {
                return false;
            }
            
            return true;
        }

        /**
         * Calculate and store the divergence of solar wind velocity for all cells in the subdomain
         */
        void SetDivSolarWindVelocity() {
            // Loop through all cells in the subdomain
            for (int CellCounter = 0; CellCounter < PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_; CellCounter++) {
                // Calculate indices for the current cell
                int nLocalNode, ii = CellCounter, i, j, k;
                
                // Get the local node index
                nLocalNode = ii / (_BLOCK_CELLS_Z_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_X_);
                ii -= nLocalNode * _BLOCK_CELLS_Z_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_X_;
                
                // Get the k index (z direction)
                k = ii / (_BLOCK_CELLS_Y_ * _BLOCK_CELLS_X_);
                ii -= k * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_X_;
                
                // Get the j index (y direction)
                j = ii / _BLOCK_CELLS_X_;
                ii -= j * _BLOCK_CELLS_X_;
                
                // Get the i index (x direction)
                i = ii;
                
                // Get the node from the block table
                cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node = PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
                auto block = node->block;
                
                // Skip if the block is NULL
                if (block == NULL) continue;
                
                // Get the cell center node
                PIC::Mesh::cDataCenterNode *cell = block->GetCenterNode(_getCenterNodeLocalNumber(i, j, k));
                
                // Skip if the cell is NULL
                if (cell == NULL) continue;
                
                // Calculate the physical coordinates of the cell center
                double x[3];
                
                // Calculate the cell center coordinates
                x[0] = node->xmin[0] + (i + 0.5) * (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_;
                x[1] = node->xmin[1] + (j + 0.5) * (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_;
                x[2] = node->xmin[2] + (k + 0.5) * (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_;
                
                // Calculate the divergence of the solar wind velocity at the cell center
                double div = GetDivSolarWindVelocity(x, node);
                
                // Save the divergence in the state vector of the cell
                *((double*) (cell->GetAssociatedDataBufferPointer() + DivSolarWindVelocityOffset)) = div;
            }
        }

        /**
         * Interpolate the divergence of solar wind velocity to a specific location
         * 
         * @param x Position vector
         * @param node Pointer to the AMR tree node
         * @return Interpolated divergence of the solar wind velocity
         */
        double InterpolateDivSolarWindVelocity(double *x, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
            double div = 0.0;
            
            // Initialize interpolation stencil
	    PIC::InterpolationRoutines::CellCentered::cStencil CenterBasedStencil;
	    PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,CenterBasedStencil);

            // Interpolate the divergence value
            for (int icell = 0; icell < CenterBasedStencil.Length; icell++) {
                // Get the data pointer and weight for this cell
                char *AssociatedDataPointer = CenterBasedStencil.cell[icell]->GetAssociatedDataBufferPointer();
                double weight = CenterBasedStencil.Weight[icell];
                
                // Get the divergence value for this cell
                double div_cell = *((double*) (AssociatedDataPointer + DivSolarWindVelocityOffset));
                
                // Add weighted contribution to interpolated value
                div += weight * div_cell;
            }
            
            return div;
        }

    } // namespace SolarWind
} // namespace SEP
