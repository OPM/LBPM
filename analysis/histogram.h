/*
 *  Generate a histogram for volumetric, interfacial and common curve properties
 *  copyright 2014, James E. McClure
 */

#define HISTOGRAM_RESOLUTION 1000

struct Histogram {
    Histogram(double v1, double v2) {
        data = new double[HISTOGRAM_RESOLUTION];
        minimum = v1;
        maximum = v2;
        delta = (maximum - minimum) / HISTOGRAM_RESOLUTION;
    }
    ~Histogram { delete *data; }
    double *data;
    double minimum, maximum, delta;

    // Adds value into the histogram
    void IncludeValue(double value, double weight) {
        idx = floor((value - min) / delta);
        if (idx > HISTOGRAM_RESOLUTION)
            ;
        else if (idx < 0)
            ;
        else
            data[idx] += weight;
    }

    // Returns the maximum value in the histogram
    void GetMax() {
        max = minimum;
        for (idx = 1; idx < HISTOGRAM_RESOLUTION; idx++) {
            if (data[idx] > max) {
                max = minimum + idx * delta;
            }
        }
    }

    // Resets the histogram to zero
    void Reset() {
        for (idx = 0; idx < HISTOGRAM_RESOLUTION; idx++)
            data[idx] = 0.0;
    }

private:
    int idx;
    double max, min;
};
