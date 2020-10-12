/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 *  Generate a histogram for volumetric, interfacial and common curve properties
 *  copyright 2014, James E. McClure
 */

#define HISTOGRAM_RESOLUTION 1000

struct Histogram{
	Histogram(double v1, double v2){
		data = new double[HISTOGRAM_RESOLUTION];
		minimum = v1;
		maximum = v2;
		delta = (maximum-minimum)/HISTOGRAM_RESOLUTION;
	}
	~Histogram{
		delete *data;
	}
	double *data;
	double minimum,maximum,delta;
	
	// Adds value into the histogram
	void IncludeValue(double value, double weight){
		idx = floor((value-min)/delta);
		if (idx > HISTOGRAM_RESOLUTION) ;
		else if (idx < 0) ;
		else data[idx] += weight;
	}
	
	// Returns the maximum value in the histogram
	void GetMax(){
		max = minimum;
		for (idx=1; idx<HISTOGRAM_RESOLUTION; idx++){
			if (data[idx] > max){
				max = minimum+idx*delta;
			}
		}
	}
	
	// Resets the histogram to zero
	void Reset(){
		for (idx=0; idx<HISTOGRAM_RESOLUTION; idx++) data[idx] = 0.0;
	}
	
private:
	int idx;
	double max,min;
};
