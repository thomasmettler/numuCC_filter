#ifndef GEOMETRYHELPER_H
#define GEOMETRYHELPER_H

#include "larcore/Geometry/Geometry.h"

// Some geometrical function
#define SMALL_NUM 0.00000001 // anything that avoids division overflow

class GeometryHelper
{
  public:
    GeometryHelper();
    ~GeometryHelper() = default;

    bool isActive(const std::vector<float> &x) const;
    bool isActive(const std::vector<double> &x) const;
    bool isActive(const float x[3]) const;
    bool isActive(const double x[3]) const;

    float dot_product_3D(const std::vector<float> &a, const std::vector<float> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
    float norm_3D(const std::vector<float> &a)
    {
        return sqrt(dot_product_3D(a, a));
    }

    float distance(const std::vector<float> &a, const std::vector<float> &b) const;
    float distance(const std::vector<double> &a, const std::vector<double> &b) const;

    // Function to calculate closest approach
    float dist3D_Line_to_Line(const std::vector<float> &l1p0, const std::vector<float> &l1p1,
                              const std::vector<float> &l2p0, const std::vector<float> &l2p1);
    // Function to calculate the cos(angle) between two lines
    float angle3D_Line_to_Line(const std::vector<float> &l1p0, const std::vector<float> &l1p1,
                               const std::vector<float> &l2p0, const std::vector<float> &l2p1);
    // Function to calculate the ratio between the summed length of the two tracks and the maximal distance between them.
    // the miximal distance should be bigger than each of the lengthts and should be comparable with the sum of the lengths if the tracks belong together.
    float fracional_maxdist(const std::vector<float> &l1p0, const std::vector<float> &l1p1,
                            const std::vector<float> &l2p0, const std::vector<float> &l2p1);

    // Function takes start and end ponts of two tracks, it returns the angle between them and the distance of the closest approach
    // Returns negative number if something went wrong.
    int brokentrack(const std::vector<float> &l1p0, const std::vector<float> &l1p1,
                    const std::vector<float> &l2p0, const std::vector<float> &l2p1,
                    float &min_dist, float &cosangle, float &total_len_fraction);

    // this function adds a track and checks if it is the broken track of a previously added one, if so it returns the index
    int add_track(const std::vector<float> &l1p0, const std::vector<float> &l1p1, const uint index);
    void clear_tracks();

  private:
    art::ServiceHandle<geo::Geometry> geo;
    // Geometry helper needs to be aware of all the tracks in the event to check for broken tracks.
    std::vector<uint> tr_index;
    std::vector<float> tr_start_x;
    std::vector<float> tr_start_y;
    std::vector<float> tr_start_z;
    std::vector<float> tr_end_x;
    std::vector<float> tr_end_y;
    std::vector<float> tr_end_z;
};

GeometryHelper::GeometryHelper()
{
    //// Check if things are set up properly:
    std::cout << std::endl;
    std::cout << "[GeometryHelper constructor] Checking set-up" << std::endl;
    std::cout << "[GeometryHelper constructor] Detector dimensions from geo: "
              << 2.0 * geo->DetHalfWidth() << ", " << geo->DetHalfHeight() << ", " << geo->DetLength() << std::endl;
}

bool GeometryHelper::isActive(const std::vector<float> &x) const
{
    if (x.size() != 3)
    {
        return false;
    }

    return this->isActive(&x[0]);
}

bool GeometryHelper::isActive(const std::vector<double> &x) const
{
    if (x.size() != 3)
    {
        return false;
    }

    return this->isActive(&x[0]);
}

bool GeometryHelper::isActive(const float x[3]) const
{
    std::vector<double> bnd = {
        0., 2.0 * geo->DetHalfWidth(),
        -1 * geo->DetHalfHeight(), geo->DetHalfHeight(),
        0., geo->DetLength()};

    bool is_x = x[0] > bnd[0] && x[0] < bnd[1];
    bool is_y = x[1] > bnd[2] && x[1] < bnd[3];
    bool is_z = x[2] > bnd[4] && x[2] < bnd[5];
    return is_x && is_y && is_z;
}

bool GeometryHelper::isActive(const double x[3]) const
{
    std::vector<double> bnd = {
        0., 2.0 * geo->DetHalfWidth(),
        -1 * geo->DetHalfHeight(), geo->DetHalfHeight(),
        0., geo->DetLength()};

    bool is_x = x[0] > bnd[0] && x[0] < bnd[1];
    bool is_y = x[1] > bnd[2] && x[1] < bnd[3];
    bool is_z = x[2] > bnd[4] && x[2] < bnd[5];
    return is_x && is_y && is_z;
}

float GeometryHelper::distance(const std::vector<float> &a,
                               const std::vector<float> &b) const
{
    if (a.size() != 3 || b.size() != 3)
    {
        return -1;
    }

    float d = 0;

    for (uint i = 0; i < 3; i++)
    {
        d += pow((a[i] - b[i]), 2);
    }

    return sqrt(d);
}

float GeometryHelper::distance(const std::vector<double> &a,
                               const std::vector<double> &b) const
{
    if (a.size() != 3 || b.size() != 3)
    {
        return -1;
    }

    float d = 0;

    for (uint i = 0; i < 3; i++)
    {
        d += pow((a[i] - b[i]), 2);
    }

    return sqrt(d);
}

// Function to calculate closest approach
float GeometryHelper::dist3D_Line_to_Line(const std::vector<float> &l1p0, const std::vector<float> &l1p1,
                                          const std::vector<float> &l2p0, const std::vector<float> &l2p1)
{
    const std::vector<float> u = {l1p1[0] - l1p0[0], l1p1[1] - l1p0[1], l1p1[2] - l1p0[2]};
    const std::vector<float> v = {l2p1[0] - l2p0[0], l2p1[1] - l2p0[1], l2p1[2] - l2p0[2]};
    const std::vector<float> w = {l1p0[0] - l2p0[0], l1p0[1] - l2p0[1], l1p0[2] - l2p0[2]};

    float a = dot_product_3D(u, u);
    float b = dot_product_3D(u, v);
    float c = dot_product_3D(v, v);
    float d = dot_product_3D(u, w);
    float e = dot_product_3D(v, w);

    float D = a * c - b * b; // D is positive
    float sc, tc;
    // compute the line parameters of the two closest points
    if (D < SMALL_NUM)
    { // the lines are almost parallel
        sc = 0.0;
        tc = (b > c ? d / b : e / c); // use the largest denominator
    }
    else
    {
        sc = (b * e - c * d) / D;
        tc = (a * e - b * d) / D;
    }

    // get the difference of the two closest points
    std::vector<float> dP = {0, 0, 0};
    for (uint i = 0; i < 3; i++)
    {
        dP[i] = w[i] + (sc * u[i]) - (tc * v[i]); // =  L1(sc) - L2(tc)
    }

    return norm_3D(dP); // return the closest distance
}

// Function to calculate the cos(angle) between two lines
float GeometryHelper::angle3D_Line_to_Line(const std::vector<float> &l1p0, const std::vector<float> &l1p1,
                                           const std::vector<float> &l2p0, const std::vector<float> &l2p1)
{
    const std::vector<float> u = {l1p1[0] - l1p0[0], l1p1[1] - l1p0[1], l1p1[2] - l1p0[2]};
    const std::vector<float> v = {l2p1[0] - l2p0[0], l2p1[1] - l2p0[1], l2p1[2] - l2p0[2]};

    return dot_product_3D(u, v) / (norm_3D(u) * norm_3D(v));
}

float GeometryHelper::fracional_maxdist(const std::vector<float> &l1p0, const std::vector<float> &l1p1,
                                        const std::vector<float> &l2p0, const std::vector<float> &l2p1)
{
    float len1 = distance(l1p0, l1p1);
    float len2 = distance(l2p0, l2p1);

    float a = distance(l2p0, l1p1);
    float b = distance(l2p0, l1p0);
    float c = distance(l2p1, l1p0);
    float d = distance(l2p1, l1p1);
    float max_dist = std::max(std::max(std::max(a, b), c), d);
    if (max_dist < len1 || max_dist < len2)
    {
        std::cout << "[GeometryHelper::fracional_maxdist] ";
        std::cout << "The maximum distance between the two segments is smaller than one of the segments itself: " << max_dist << std::endl;
        return 0;
    }
    else
    {
        return max_dist / (len1 + len2);
    }
}

int GeometryHelper::brokentrack(const std::vector<float> &l1p0, const std::vector<float> &l1p1,
                                const std::vector<float> &l2p0, const std::vector<float> &l2p1,
                                float &min_dist, float &cosangle, float &total_len_fraction)
{
    //parameters to decide if tracks are actually broken brother:
    float mincos = 0.9;
    float maxdist = 0.1;
    float frac_min = 0.75;
    float frac_max = 1.25;

    cosangle = angle3D_Line_to_Line(l1p0, l1p1,
                                    l2p0, l2p1);
    min_dist = dist3D_Line_to_Line(l1p0, l1p1,
                                   l2p0, l2p1);

    total_len_fraction = fracional_maxdist(l1p0, l1p1,
                                           l2p0, l2p1);

    std::cout << "[brokentrack] "
              << "track 1 start: (" << int(l1p0[0]) << ", " << int(l1p0[1]) << ", " << int(l1p0[2]) << ")";
    std::cout << " end: (" << int(l1p1[0]) << ", " << int(l1p1[1]) << ", " << int(l1p1[2]) << ")" << std::endl;
    std::cout << "[brokentrack] "
              << "track 2 start: (" << int(l2p0[0]) << ", " << int(l2p0[1]) << ", " << int(l2p0[2]) << ")";
    std::cout << " end: (" << int(l2p1[0]) << ", " << int(l2p1[1]) << ", " << int(l2p1[2]) << ")" << std::endl;

    std::cout << "[brokentrack] "
              << "(cosangle, min_dist, length_fraction): " << cosangle << ", " << int(l2p1[1]) << ", " << int(l2p1[2]) << ")" << std::endl;

    if (std::abs(cosangle) > mincos && min_dist < maxdist && total_len_fraction > frac_min && total_len_fraction < frac_max)
    {
        std::cout << "[brokentrack] pair of broken tracks identified!" << std::endl;
        return 0;
    }
    else
    {
        return -1;
    }
}

void GeometryHelper::clear_tracks()
{
    tr_index.clear();
    tr_start_x.clear();
    tr_start_y.clear();
    tr_start_z.clear();
    tr_end_x.clear();
    tr_end_y.clear();
    tr_end_z.clear();
}

// this function adds a track and checks if it is the broken track of a previously added one, if so it returns the index
int GeometryHelper::add_track(const std::vector<float> &l1p0, const std::vector<float> &l1p1, const uint index)
{
    int result = -1;
    // Does it have brothers?
    uint n_tr = tr_index.size();
    for (uint i = 0; i < n_tr; ++i)
    {
        float min_dist;
        float cosangle;
        float total_len_fraction;

        const std::vector<float> l2p0 = {tr_start_x[i], tr_start_y[i], tr_start_z[i]};
        const std::vector<float> l2p1 = {tr_end_x[i], tr_end_y[i], tr_end_z[i]};
        bool is_broken = brokentrack(l1p0, l1p1, l2p0, l2p1, min_dist, cosangle, total_len_fraction);
        if (is_broken)
        {
            result = tr_index[i];
            std::cout << "[add_track] Track " << index << " is matched with track " << i << std::endl;
        }
    }

    // Add the track
    tr_index.push_back(index);
    tr_start_x.push_back(l1p0[0]);
    tr_start_y.push_back(l1p0[1]);
    tr_start_z.push_back(l1p0[2]);
    tr_end_x.push_back(l1p1[0]);
    tr_end_y.push_back(l1p1[1]);
    tr_end_z.push_back(l1p1[2]);
    return result;
}

#endif