#ifndef SEP_SOLAR_WIND_H
#define SEP_SOLAR_WIND_H

/**
 * @file SolarWind.h
 * @brief Header file defining solar wind related functions
 */

namespace SEP {
    namespace SolarWind {
        /**
         * Offset for solar wind velocity divergence in the data buffer
         */
        extern int DivSolarWindVelocityOffset;

        /**
         * Calculates the solar wind velocity at a given position.
         * 
         * @param v Output array for velocity components
         * @param x Position vector
         * @return The magnitude of the solar wind velocity
         */
        double GetSolarWindVelocity(double *v, double *x);
	void Init();

        /**
         * Calculate and store the divergence of solar wind velocity for all cells in the subdomain
         */
        void SetDivSolarWindVelocity();  

        /**
         * Interpolate the divergence of solar wind velocity to a specific location
         * 
         * @param x Position vector
         * @param node Pointer to the AMR tree node
         * @return Interpolated divergence of the solar wind velocity
         */
        double InterpolateDivSolarWindVelocity(double *x, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node); 

        /**
         * Check if a point is valid (not outside domain or inside the Sun)
         * 
         * @param x Position vector to check
         * @return true if valid, false otherwise
         */
        bool IsPointValid(double *x);
        
    } // namespace SolarWind
} // namespace SEP

#endif // SEP_SOLAR_WIND_H
