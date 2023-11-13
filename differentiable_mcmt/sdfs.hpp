#include <cmath>


namespace SDF{
    double sphere_sdf(double x, double y, double z)
    {
        // Center of the sphere
        double centerX = 0.5;
        double centerY = 0.5;
        double centerZ = 0.5;

        centerX = 0.0;
        centerY = 0.0;
        centerZ = 0.0;


        // Radius of the sphere
        double radius = 0.3; // Adjust the radius as needed (should be less than 0.5 to fit within [0, 1])

        // Calculate the distance from the point to the center of the sphere
        double distance = std::sqrt((x - centerX) * (x - centerX) + (y - centerY) * (y - centerY) + (z - centerZ) * (z - centerZ));

        // Return the signed distance (negative inside the sphere, positive outside the sphere)
        return distance - radius;
    }

    double sdBox(double px, double py, double pz)
    {

        px -= 0.5;
        py -= 0.5;
        pz -= 0.5;

        double bx = 0.3;
        double by = 0.3;
        double bz = 0.3;

        double qx = abs(px) - bx;
        double qy = abs(py) - by;
        double qz = abs(pz) - bz;

        double max_q = std::max(std::max(qx, qy), qz);
        double length_max_q = std::sqrt(std::max(qx, 0.0)*std::max(qx, 0.0) + std::max(qy, 0.0)*std::max(qy, 0.0) + std::max(qz, 0.0)*std::max(qz, 0.0));

        return length_max_q + std::min(max_q, 0.0);
    }
    
}