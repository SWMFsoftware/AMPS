// WeightedParticleMergingSplitting.cpp
// This program implements Weighted Particle Merging and Splitting algorithms tailored for particles characterized by
// one spatial coordinate (s) and two velocity coordinates (vParallel and vNormal).
// The algorithms conserve momentum, report energy discrepancies,
// and prioritize merging/splitting based on particle weights and bin population.

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <array>
#include <algorithm>
#include <iomanip>

// Placeholder namespaces and functions for PB and PIC
// In actual implementation, these should interface with the simulation's particle data structures

namespace PB {
    // Retrieve the spatial coordinate 's' of particle p
    double GetFieldLineCoord(long int p) {
        // Placeholder implementation
        // Replace with actual data retrieval
        // Example:
        // return particleData[p].s;
        return static_cast<double>(p % 1000) / 1000.0; // Example: s in [0,1)
    }

    // Retrieve the parallel velocity component of particle p
    double GetVParallel(long int p) {
        // Placeholder implementation
        // Replace with actual data retrieval
        // Example:
        // return particleData[p].vParallel;
        return 10.0 + 0.1 * p; // Example velocities
    }

    // Retrieve the normal velocity component of particle p
    double GetVNormal(long int p) {
        // Placeholder implementation
        // Replace with actual data retrieval
        // Example:
        // return particleData[p].vNormal;
        return 5.0 + 0.05 * p; // Example velocities
    }

    // Retrieve the weight (mass) of particle p
    double GetWeight(long int p) {
        // Placeholder implementation
        // Replace with actual data retrieval
        // Example:
        // return particleData[p].w;
        return 1.0 + 0.05 * p; // Example weight
    }

    // Remove particle p from the linked list
    void RemoveParticle(long int p) {
        // Placeholder implementation
        // Implement actual particle removal from the linked list
        // Example:
        // linkedList.remove(p);
        // For demonstration, print removal
        std::cout << "Removing particle ID " << p << "\n";
    }

    // Add a new particle to the linked list and return its identifier
    long int AddParticle(double s, double vParallel, double vNormal, double w) {
        // Placeholder implementation
        // Implement actual particle addition and return new particle identifier
        // Example:
        // linkedList.add(newParticle);
        // return newParticle.id;
        static long int new_id = 10000; // Starting ID for new particles
        new_id++;
        // For demonstration, print the new particle's properties
        std::cout << "Adding merged/split particle ID " << new_id << " with properties:\n";
        std::cout << "  s: " << s << "\n";
        std::cout << "  vParallel: " << vParallel << " m/s\n";
        std::cout << "  vNormal: " << vNormal << " m/s\n";
        std::cout << "  Weight: " << w << " kg\n\n";
        return new_id;
    }
}

namespace PIC {
    namespace ParticleBuffer {
        // Retrieve the next particle in the linked list
        long int GetNext(long int p) {
            // Placeholder implementation
            // Replace with actual method to get the next particle
            // For demonstration, return -1 to indicate end of list
            return -1;
        }
    }
}

// Structure to hold particle data for binning
struct ParticleData {
    double s;          // Spatial coordinate
    double vParallel;  // Parallel velocity component
    double vNormal;    // Normal velocity component
    double w;          // Weight (mass)
    long int id;       // Particle identifier
};

// Function to compute a unique key for a bin based on its spatial and velocity indices
long long computeBinKey(int is, int ivParallel, int ivNormal, 
                        int numBinsSpatial, int numBinsVParallel, int numBinsVNormal) {
    // Using prime number multipliers to generate a unique key
    const long long prime1 = 73856093;
    const long long prime2 = 19349663;
    const long long prime3 = 83492791;

    return static_cast<long long>(is) * prime1 ^ 
           static_cast<long long>(ivParallel) * prime2 ^ 
           static_cast<long long>(ivNormal) * prime3;
}

/// \brief Performs Weighted Particle Merging on a linked list of particles.
///
/// This function groups particles into spatial and velocity bins and merges them within each combined bin
/// based on their weights. The merging process prioritizes combining particles with lower weights first, 
/// two at a time, until the total number of particles is within the specified range [nParticleRangeMin, nParticleRangeMax].
/// Momentum is conserved exactly, and energy discrepancies are reported.
///
/// \param head The identifier of the first particle in the linked list.
/// \param numBinsSpatial The number of spatial bins along the 's' axis.
/// \param numBinsVParallel The number of velocity bins along the 'vParallel' axis.
/// \param numBinsVNormal The number of velocity bins along the 'vNormal' axis.
/// \param sRange The range of spatial coordinate 's' (typically [0,1)).
/// \param vParallelRange Reference to a variable that will store the dynamically calculated maximum absolute vParallel.
/// \param vNormalRange Reference to a variable that will store the dynamically calculated maximum absolute vNormal.
/// \param nParticleRangeMin The minimum desired number of particles after merging.
/// \param nParticleRangeMax The maximum desired number of particles after merging.
/// \param mergeThreshold Minimum number of particles in a bin to perform merging (default is 2).
void WeightedParticleMerging(long int head, int numBinsSpatial, int numBinsVParallel, int numBinsVNormal,
                             double sRange, double &vParallelRange, double &vNormalRange,
                             int nParticleRangeMin, int nParticleRangeMax,
                             int mergeThreshold = 2) {
    // Initial bin sizes (will be updated after dynamic V_MAX calculation)
    double binSizeSpatial = sRange / numBinsSpatial;
    double binSizeVParallel = 0.0;
    double binSizeVNormal = 0.0;

    // Data structure to hold particles in each combined spatial-velocity bin
    // Key: combined bin key, Value: vector of ParticleData
    std::unordered_map<long long, std::vector<ParticleData>> bins;

    // Variables to determine V_MAX dynamically
    double dynamicVParallelMax = 0.0;
    double dynamicVNormalMax = 0.0;

    // First pass: Traverse all particles to determine dynamic V_MAX
    long int p_temp = head;
    int totalParticles = 0;
    while (p_temp != -1) {
        double vParallel = PB::GetVParallel(p_temp);
        double vNormal = PB::GetVNormal(p_temp);

        dynamicVParallelMax = std::max(dynamicVParallelMax, std::abs(vParallel));
        dynamicVNormalMax = std::max(dynamicVNormalMax, std::abs(vNormal));

        totalParticles++;
        p_temp = PIC::ParticleBuffer::GetNext(p_temp);
    }

    // Update velocity ranges based on dynamic V_MAX
    vParallelRange = dynamicVParallelMax;
    vNormalRange = dynamicVNormalMax;

    // Calculate bin sizes based on dynamic V_MAX
    binSizeVParallel = (2.0 * vParallelRange) / numBinsVParallel; // Assuming vParallel ∈ [-V_MAX, +V_MAX]
    binSizeVNormal = (2.0 * vNormalRange) / numBinsVNormal;       // Assuming vNormal ∈ [-V_MAX, +V_MAX]

    // Reset particle traversal
    p_temp = head;

    // Traverse the linked list and assign particles to combined spatial-velocity bins
    while (p_temp != -1) {
        // Retrieve particle properties
        double s = PB::GetFieldLineCoord(p_temp);
        double vParallel = PB::GetVParallel(p_temp);
        double vNormal = PB::GetVNormal(p_temp);
        double w = PB::GetWeight(p_temp);

        // Determine spatial bin index based on 's'
        int is = static_cast<int>(std::floor(s / binSizeSpatial));
        // Clamp the index to [0, numBinsSpatial - 1]
        if (is >= numBinsSpatial) is = numBinsSpatial - 1;
        if (is < 0) is = 0;

        // Determine velocity bin indices based on 'vParallel' and 'vNormal'
        int ivParallel = static_cast<int>(std::floor((vParallel + vParallelRange) / binSizeVParallel));
        int ivNormal = static_cast<int>(std::floor((vNormal + vNormalRange) / binSizeVNormal));

        // Clamp the indices to valid ranges
        if (ivParallel >= numBinsVParallel) ivParallel = numBinsVParallel - 1;
        if (ivParallel < 0) ivParallel = 0;
        if (ivNormal >= numBinsVNormal) ivNormal = numBinsVNormal - 1;
        if (ivNormal < 0) ivNormal = 0;

        // Compute combined bin key
        long long binKey = computeBinKey(is, ivParallel, ivNormal, numBinsSpatial, numBinsVParallel, numBinsVNormal);

        // Store particle data in the corresponding bin
        ParticleData pdata;
        pdata.s = s;
        pdata.vParallel = vParallel;
        pdata.vNormal = vNormal;
        pdata.w = w;
        pdata.id = p_temp;

        bins[binKey].push_back(pdata);

        // Move to the next particle
        p_temp = PIC::ParticleBuffer::GetNext(p_temp);
    }

    // Check if merging is needed
    if (totalParticles <= nParticleRangeMax && totalParticles >= nParticleRangeMin) {
        std::cout << "No merging needed. Total particles (" << totalParticles << ") are within the desired range [" 
                  << nParticleRangeMin << ", " << nParticleRangeMax << "].\n";
    }

    // If totalParticles < nParticleRangeMin, merging might not be necessary or desired.
    // Depending on simulation needs, implement splitting here. For now, we focus on reducing particles.

    // Calculate how many particles need to be removed
    int particlesToRemove = (totalParticles > nParticleRangeMax) ? (totalParticles - nParticleRangeMax) : 0;

    if (particlesToRemove > 0) {
        // Collect all bins that can be merged (i.e., have at least 'mergeThreshold' particles)
        std::vector<std::pair<long long, std::vector<ParticleData>>> mergableBins;
        for (const auto &bin : bins) {
            if (bin.second.size() >= mergeThreshold) {
                mergableBins.emplace_back(bin);
            }
        }

        // Sort the mergable bins based on the average weight ascending (prioritize lower-weight bins)
        std::sort(mergableBins.begin(), mergableBins.end(),
                  [&](const std::pair<long long, std::vector<ParticleData>> &a,
                      const std::pair<long long, std::vector<ParticleData>> &b) -> bool {
                        double sumWeightA = 0.0;
                        for (const auto &p : a.second) sumWeightA += p.w;
                        double avgWeightA = sumWeightA / a.second.size();

                        double sumWeightB = 0.0;
                        for (const auto &p : b.second) sumWeightB += p.w;
                        double avgWeightB = sumWeightB / b.second.size();

                        return avgWeightA < avgWeightB;
                  });

        // Variables to track total energy before and after merging
        double totalEnergyBefore = 0.0;
        double totalEnergyAfter = 0.0;

        // Start merging process
        for (auto &bin : mergableBins) {
            if (particlesToRemove <= 0) break; // Desired range achieved

            auto &particles = bin.second;

            // Sort particles in the bin by weight ascending
            std::sort(particles.begin(), particles.end(),
                      [&](const ParticleData &a, const ParticleData &b) -> bool {
                          return a.w < b.w;
                      });

            // Continue merging within this bin as long as:
            // - There are enough particles to merge
            // - Merging reduces the required number of particles
            while (particles.size() >= 2 && particlesToRemove > 0) {
                // Select the two lowest-weight particles
                ParticleData p1 = particles[0];
                ParticleData p2 = particles[1];

                // Calculate total weight and momentum
                double mergedWeight = p1.w + p2.w;
                double totalMomentumParallel = (p1.w * p1.vParallel) + (p2.w * p2.vParallel);
                double totalMomentumNormal = (p1.w * p1.vNormal) + (p2.w * p2.vNormal);
                double mergedSWeighted = (p1.s * p1.w) + (p2.s * p2.w);
                double mergedEnergyBefore = 0.5 * p1.w * (p1.vParallel * p1.vParallel + p1.vNormal * p1.vNormal) +
                                            0.5 * p2.w * (p2.vParallel * p2.vParallel + p2.vNormal * p2.vNormal);

                // Compute merged particle's velocity to conserve momentum
                double mergedVParallel = totalMomentumParallel / mergedWeight;
                double mergedVNormal = totalMomentumNormal / mergedWeight;

                // Compute merged particle's kinetic energy
                double mergedEnergyAfter = 0.5 * mergedWeight * (mergedVParallel * mergedVParallel + mergedVNormal * mergedVNormal);

                // Accumulate energies for reporting
                totalEnergyBefore += mergedEnergyBefore;
                totalEnergyAfter += mergedEnergyAfter;

                // Compute weighted average spatial coordinate 's'
                double mergedS = mergedSWeighted / mergedWeight;

                // Create a new merged particle
                long int mergedParticle = PB::AddParticle(mergedS, mergedVParallel, mergedVNormal, mergedWeight);

                // Remove original particles from the linked list
                PB::RemoveParticle(p1.id);
                PB::RemoveParticle(p2.id);
                totalParticles -= 2;
                particlesToRemove -= 2;

                // Remove the first two particles from the bin's particle list
                particles.erase(particles.begin(), particles.begin() + 2);
            }
        }

        // After merging, check if the desired range is achieved
        if (particlesToRemove > 0) {
            std::cout << "Warning: Unable to reduce particles to the desired maximum range.\n";
            std::cout << "Remaining particles to remove: " << particlesToRemove << "\n";
        }

        // Report energy discrepancies
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Total Kinetic Energy Before Merging: " << totalEnergyBefore << " J\n";
        std::cout << "Total Kinetic Energy After Merging:  " << totalEnergyAfter << " J\n";
        std::cout << "Energy Discrepancy:                 " << (mergedEnergyAfter - mergedEnergyBefore) << " J\n\n";

        std::cout << "Weighted Particle Merging completed.\n";
    }

/// \brief Performs Weighted Particle Splitting on a linked list of particles.
///
/// This function groups particles into spatial and velocity bins and splits them within each combined bin
/// based on their weights. The splitting process prioritizes splitting particles with larger weights first 
/// in the most populated bins, until the total number of particles is within the specified range [nParticleRangeMin, nParticleRangeMax].
/// Momentum is conserved exactly, and energy discrepancies are reported.
///
/// \param head The identifier of the first particle in the linked list.
/// \param numBinsSpatial The number of spatial bins along the 's' axis.
/// \param numBinsVParallel The number of velocity bins along the 'vParallel' axis.
/// \param numBinsVNormal The number of velocity bins along the 'vNormal' axis.
/// \param sRange The range of spatial coordinate 's' (typically [0,1)).
/// \param vParallelRange Reference to a variable that will store the dynamically calculated maximum absolute vParallel.
/// \param vNormalRange Reference to a variable that will store the dynamically calculated maximum absolute vNormal.
/// \param nParticleRangeMin The minimum desired number of particles after splitting.
/// \param nParticleRangeMax The maximum desired number of particles after splitting.
/// \param splitThreshold Minimum number of particles in a bin to perform splitting (default is 2).
void WeightedParticleSplitting(long int head, int numBinsSpatial, int numBinsVParallel, int numBinsVNormal,
                               double sRange, double &vParallelRange, double &vNormalRange,
                               int nParticleRangeMin, int nParticleRangeMax,
                               int splitThreshold = 2) {
    // Initial bin sizes (will be updated after dynamic V_MAX calculation)
    double binSizeSpatial = sRange / numBinsSpatial;
    double binSizeVParallel = 0.0;
    double binSizeVNormal = 0.0;

    // Data structure to hold particles in each combined spatial-velocity bin
    // Key: combined bin key, Value: vector of ParticleData
    std::unordered_map<long long, std::vector<ParticleData>> bins;

    // Variables to determine V_MAX dynamically
    double dynamicVParallelMax = 0.0;
    double dynamicVNormalMax = 0.0;

    // First pass: Traverse all particles to determine dynamic V_MAX
    long int p_temp = head;
    int totalParticles = 0;
    while (p_temp != -1) {
        double vParallel = PB::GetVParallel(p_temp);
        double vNormal = PB::GetVNormal(p_temp);

        dynamicVParallelMax = std::max(dynamicVParallelMax, std::abs(vParallel));
        dynamicVNormalMax = std::max(dynamicVNormalMax, std::abs(vNormal));

        totalParticles++;
        p_temp = PIC::ParticleBuffer::GetNext(p_temp);
    }

    // Update velocity ranges based on dynamic V_MAX
    vParallelRange = dynamicVParallelMax;
    vNormalRange = dynamicVNormalMax;

    // Calculate bin sizes based on dynamic V_MAX
    binSizeVParallel = (2.0 * vParallelRange) / numBinsVParallel; // Assuming vParallel ∈ [-V_MAX, +V_MAX]
    binSizeVNormal = (2.0 * vNormalRange) / numBinsVNormal;       // Assuming vNormal ∈ [-V_MAX, +V_MAX]

    // Reset particle traversal
    p_temp = head;

    // Traverse the linked list and assign particles to combined spatial-velocity bins
    while (p_temp != -1) {
        // Retrieve particle properties
        double s = PB::GetFieldLineCoord(p_temp);
        double vParallel = PB::GetVParallel(p_temp);
        double vNormal = PB::GetVNormal(p_temp);
        double w = PB::GetWeight(p_temp);

        // Determine spatial bin index based on 's'
        int is = static_cast<int>(std::floor(s / binSizeSpatial));
        // Clamp the index to [0, numBinsSpatial - 1]
        if (is >= numBinsSpatial) is = numBinsSpatial - 1;
        if (is < 0) is = 0;

        // Determine velocity bin indices based on 'vParallel' and 'vNormal'
        int ivParallel = static_cast<int>(std::floor((vParallel + vParallelRange) / binSizeVParallel));
        int ivNormal = static_cast<int>(std::floor((vNormal + vNormalRange) / binSizeVNormal));

        // Clamp the indices to valid ranges
        if (ivParallel >= numBinsVParallel) ivParallel = numBinsVParallel - 1;
        if (ivParallel < 0) ivParallel = 0;
        if (ivNormal >= numBinsVNormal) ivNormal = numBinsVNormal - 1;
        if (ivNormal < 0) ivNormal = 0;

        // Compute combined bin key
        long long binKey = computeBinKey(is, ivParallel, ivNormal, numBinsSpatial, numBinsVParallel, numBinsVNormal);

        // Store particle data in the corresponding bin
        ParticleData pdata;
        pdata.s = s;
        pdata.vParallel = vParallel;
        pdata.vNormal = vNormal;
        pdata.w = w;
        pdata.id = p_temp;

        bins[binKey].push_back(pdata);

        // Move to the next particle
        p_temp = PIC::ParticleBuffer::GetNext(p_temp);
    }

    // Check if splitting is needed
    if (totalParticles <= nParticleRangeMax && totalParticles >= nParticleRangeMin) {
        std::cout << "No splitting needed. Total particles (" << totalParticles << ") are within the desired range [" 
                  << nParticleRangeMin << ", " << nParticleRangeMax << "].\n";
    }

    // Calculate how many particles need to be added
    int particlesToAdd = (totalParticles < nParticleRangeMin) ? (nParticleRangeMin - totalParticles) : 0;

    if (particlesToAdd > 0) {
        // Collect all bins that can be split (i.e., have at least 'splitThreshold' particles)
        std::vector<std::pair<long long, std::vector<ParticleData>>> splittableBins;
        for (const auto &bin : bins) {
            if (bin.second.size() >= splitThreshold) {
                splittableBins.emplace_back(bin);
            }
        }

        // Sort the splittable bins based on the number of particles descending (prioritize more populated bins)
        std::sort(splittableBins.begin(), splittableBins.end(),
                  [&](const std::pair<long long, std::vector<ParticleData>> &a,
                      const std::pair<long long, std::vector<ParticleData>> &b) -> bool {
                        return a.second.size() > b.second.size();
                  });

        // Variables to track total energy before and after splitting
        double totalEnergyBefore = 0.0;
        double totalEnergyAfter = 0.0;

        // Start splitting process
        for (auto &bin : splittableBins) {
            if (particlesToAdd <= 0) break; // Desired range achieved

            auto &particles = bin.second;

            // Sort particles in the bin by weight descending (prioritize heavier particles)
            std::sort(particles.begin(), particles.end(),
                      [&](const ParticleData &a, const ParticleData &b) -> bool {
                          return a.w > b.w;
                      });

            // Continue splitting within this bin as long as:
            // - There are enough particles to split
            // - Splitting reduces the required number of particles
            while (particles.size() >= splitThreshold && particlesToAdd > 0) {
                // Select the highest-weight particle
                ParticleData p_to_split = particles[0];

                // Calculate properties for the two new particles
                double newWeight = p_to_split.w / 2.0;

                // Ensure that the weight is above a minimum threshold to avoid infinite splitting
                if (newWeight < 1e-6) { // Example threshold
                    std::cout << "Warning: Particle ID " << p_to_split.id << " has too small weight to split.\n";
                    break;
                }

                // Compute kinetic energy before splitting
                double energyBefore = 0.5 * p_to_split.w * (p_to_split.vParallel * p_to_split.vParallel + p_to_split.vNormal * p_to_split.vNormal);

                // Create two new particles with half the weight each
                long int newParticle1 = PB::AddParticle(p_to_split.s, p_to_split.vParallel, p_to_split.vNormal, newWeight);
                long int newParticle2 = PB::AddParticle(p_to_split.s, p_to_split.vParallel, p_to_split.vNormal, newWeight);

                // Compute kinetic energy after splitting
                double energyAfter = 0.5 * newWeight * (p_to_split.vParallel * p_to_split.vParallel + p_to_split.vNormal * p_to_split.vNormal) * 2;

                // Accumulate energies for reporting
                totalEnergyBefore += energyBefore;
                totalEnergyAfter += energyAfter;

                // Remove the original particle from the linked list
                PB::RemoveParticle(p_to_split.id);
                totalParticles--;
                particlesToAdd -= 2;

                // Remove the split particle from the bin's particle list
                particles.erase(particles.begin());

                // Check if we've added enough particles
                if (particlesToAdd <= 0) break;
            }
        }

        // After splitting, check if the desired range is achieved
        if (particlesToAdd > 0) {
            std::cout << "Warning: Unable to increase particles to the desired minimum range.\n";
            std::cout << "Remaining particles to add: " << particlesToAdd << "\n";
        }

        // Report energy discrepancies
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Total Kinetic Energy Before Splitting: " << totalEnergyBefore << " J\n";
        std::cout << "Total Kinetic Energy After Splitting:  " << totalEnergyAfter << " J\n";
        std::cout << "Energy Discrepancy:                   " << (totalEnergyAfter - totalEnergyBefore) << " J\n\n";

        std::cout << "Weighted Particle Splitting completed.\n";
    }

    /*EXAMPLE:
/// \brief Example usage of the WeightedParticleMerging and WeightedParticleSplitting functions.
int main() {
    // Define binning parameters
    int numBinsSpatial = 10;       // Number of spatial bins along 's'
    int numBinsVParallel = 10;     // Number of velocity bins along 'vParallel'
    int numBinsVNormal = 10;       // Number of velocity bins along 'vNormal'
    double sRange = 1.0;           // Spatial coordinate range [0,1)

    // Initial velocity ranges (will be updated dynamically)
    double vParallelRange = 100.0; // Initial guess; will be updated based on particle velocities
    double vNormalRange = 100.0;   // Initial guess; will be updated based on particle velocities

    // Define desired particle range
    int nParticleRangeMin = 50;    // Minimum desired number of particles
    int nParticleRangeMax = 100;   // Maximum desired number of particles

    // Identifier of the first particle in the linked list
    long int head = 1; // Starting with particle ID 1

    // Perform Weighted Particle Merging with Prioritized Pairwise Low-Weight Merging
    WeightedParticleMerging(head, numBinsSpatial, numBinsVParallel, numBinsVNormal, 
                            sRange, vParallelRange, vNormalRange, 
                            nParticleRangeMin, nParticleRangeMax, 2); // Merge bins with 2 or more particles

    // Perform Weighted Particle Splitting with Prioritized Pairwise High-Weight Splitting
    WeightedParticleSplitting(head, numBinsSpatial, numBinsVParallel, numBinsVNormal, 
                              sRange, vParallelRange, vNormalRange, 
                              nParticleRangeMin, nParticleRangeMax, 2); // Split bins with 2 or more particles

    return 0;
}
*/ 
