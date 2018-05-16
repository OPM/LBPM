#include "analysis/eikonal.h"
#include "analysis/imfilter.h"




static inline float minmod(float &a, float &b)
{
    float value = a;
    if     ( a*b < 0.0)
        value=0.0;
    else if (fabs(a) > fabs(b))
        value = b;
    return value;
}


static inline double minmod(double &a, double &b){

	double value;

	value = a;
	if 	( a*b < 0.0)	    value=0.0;
	else if (fabs(a) > fabs(b)) value = b;

	return value;
}


template<class TYPE>
static inline MPI_Datatype getType( );
template<> inline MPI_Datatype getType<float>() { return MPI_FLOAT; }
template<> inline MPI_Datatype getType<double>() { return MPI_DOUBLE; }


/******************************************************************
* Solve the eikonal equation                                      *
******************************************************************/
template<class TYPE>
TYPE Eikonal( Array<TYPE> &Distance, const Array<char> &ID, const Domain &Dm, const int timesteps)
{

    /*
    * This routine converts the data in the Distance array to a signed distance
    * by solving the equation df/dt = sign(1-|grad f|), where Distance provides
    * the values of f on the mesh associated with domain Dm
    * It has been tested with segmented data initialized to values [-1,1]
    * and will converge toward the signed distance to the surface bounding the associated phases
    *
    * Reference:
    * Min C (2010) On reinitializing level set functions, Journal of Computational Physics 229
    */

    int xdim = Dm.Nx-2;
    int ydim = Dm.Ny-2;
    int zdim = Dm.Nz-2;
    fillHalo<TYPE> fillData(Dm.Comm, Dm.rank_info,xdim,ydim,zdim,1,1,1,0,1);

    // Arrays to store the second derivatives
    Array<TYPE> Dxx(Dm.Nx,Dm.Ny,Dm.Nz);
    Array<TYPE> Dyy(Dm.Nx,Dm.Ny,Dm.Nz);
    Array<TYPE> Dzz(Dm.Nx,Dm.Ny,Dm.Nz);
    Array<int8_t> sign(ID.size());
    for (size_t i=0; i<sign.length(); i++)
        sign(i) = ID(i) == 0 ? -1:1;

    int count = 0;
    double dt = 0.1;
    double GlobalVar = 0.0;
    while (count < timesteps){

        // Communicate the halo of values
        fillData.fill(Distance);

        // Compute second order derivatives
        for (int k=1;k<Dm.Nz-1;k++){
            for (int j=1;j<Dm.Ny-1;j++){
                for (int i=1;i<Dm.Nx-1;i++){
                    Dxx(i,j,k) = Distance(i+1,j,k) + Distance(i-1,j,k) - 2*Distance(i,j,k);
                    Dyy(i,j,k) = Distance(i,j+1,k) + Distance(i,j-1,k) - 2*Distance(i,j,k);
                    Dzz(i,j,k) = Distance(i,j,k+1) + Distance(i,j,k-1) - 2*Distance(i,j,k);
                }
            }
        }
        fillData.fill(Dxx);
        fillData.fill(Dyy);
        fillData.fill(Dzz);

        double LocalMax = 0.0;
        double LocalVar = 0.0;
        // Execute the next timestep
        for (int k=1;k<Dm.Nz-1;k++){
            for (int j=1;j<Dm.Ny-1;j++){
                for (int i=1;i<Dm.Nx-1;i++){
                    double s = sign(i,j,k);

                    // local second derivative terms
                    double Dxxp = minmod(Dxx(i,j,k),Dxx(i+1,j,k));
                    double Dyyp = minmod(Dyy(i,j,k),Dyy(i,j+1,k));
                    double Dzzp = minmod(Dzz(i,j,k),Dzz(i,j,k+1));
                    double Dxxm = minmod(Dxx(i,j,k),Dxx(i-1,j,k));
                    double Dyym = minmod(Dyy(i,j,k),Dyy(i,j-1,k));
                    double Dzzm = minmod(Dzz(i,j,k),Dzz(i,j,k-1));

                    /* //............Compute upwind derivatives ...................
                    Dxp = Distance(i+1,j,k) - Distance(i,j,k) + 0.5*Dxxp;
                    Dyp = Distance(i,j+1,k) - Distance(i,j,k) + 0.5*Dyyp;
                    Dzp = Distance(i,j,k+1) - Distance(i,j,k) + 0.5*Dzzp;
                    Dxm = Distance(i,j,k) - Distance(i-1,j,k) + 0.5*Dxxm;
                    Dym = Distance(i,j,k) - Distance(i,j-1,k) + 0.5*Dyym;
                    Dzm = Distance(i,j,k) - Distance(i,j,k-1) + 0.5*Dzzm;
                    */
                    double Dxp = Distance(i+1,j,k)- Distance(i,j,k) - 0.5*Dxxp;
                    double Dyp = Distance(i,j+1,k)- Distance(i,j,k) - 0.5*Dyyp;
                    double Dzp = Distance(i,j,k+1)- Distance(i,j,k) - 0.5*Dzzp;
                    double Dxm = Distance(i,j,k) - Distance(i-1,j,k) + 0.5*Dxxm;
                    double Dym = Distance(i,j,k) - Distance(i,j-1,k) + 0.5*Dyym;
                    double Dzm = Distance(i,j,k) - Distance(i,j,k-1) + 0.5*Dzzm;
                    // Compute upwind derivatives for Godunov Hamiltonian
                    double Dx, Dy, Dz;
                    if ( s < 0.0){
                        if (Dxp + Dxm > 0.f)
                            Dx = Dxp*Dxp;
                        else
                            Dx = Dxm*Dxm;
                        if (Dyp + Dym > 0.f)
                            Dy = Dyp*Dyp;
                        else
                            Dy = Dym*Dym;
                        if (Dzp + Dzm > 0.f)
                            Dz = Dzp*Dzp;
                        else
                            Dz = Dzm*Dzm;
                    }
                    else{
                        if (Dxp + Dxm < 0.f)
                            Dx = Dxp*Dxp;
                        else
                            Dx = Dxm*Dxm;
                        if (Dyp + Dym < 0.f)
                            Dy = Dyp*Dyp;
                        else
                            Dy = Dym*Dym;
                        if (Dzp + Dzm < 0.f)
                            Dz = Dzp*Dzp;
                        else
                            Dz = Dzm*Dzm;
                    }

                    //Dx = max(Dxp*Dxp,Dxm*Dxm);
                    //Dy = max(Dyp*Dyp,Dym*Dym);
                    //Dz = max(Dzp*Dzp,Dzm*Dzm);

                    double norm = sqrt(Dx + Dy + Dz);
                    if (norm > 1.0)
                        norm = 1.0;

                    Distance(i,j,k) += dt*s*(1.0 - norm);
                    LocalVar +=  dt*s*(1.0 - norm);

                    if (fabs(dt*s*(1.0 - norm)) > LocalMax)
                        LocalMax = fabs(dt*s*(1.0 - norm));
                }
            }
        }

        double GlobalMax;
        MPI_Allreduce(&LocalVar,&GlobalVar,1,MPI_DOUBLE,MPI_SUM,Dm.Comm);
        MPI_Allreduce(&LocalMax,&GlobalMax,1,MPI_DOUBLE,MPI_MAX,Dm.Comm);
        GlobalVar /= (Dm.Nx-2)*(Dm.Ny-2)*(Dm.Nz-2)*Dm.nprocx()*Dm.nprocy()*Dm.nprocz();
        count++;


        if (count%50 == 0 && Dm.rank()==0 ){
            printf("Time=%i, Max variation=%f, Global variation=%f \n",count,GlobalMax,GlobalVar);
            fflush(stdout);
        }

        if (fabs(GlobalMax) < 1e-5){
            if (Dm.rank()==0)
                printf("Exiting with max tolerance of 1e-5 \n");
            count=timesteps;
        }
    }
    return GlobalVar;
}
template float Eikonal<float>( Array<float>&, const Array<char>&, const Domain&, int );
template double Eikonal<double>( Array<double>&, const Array<char>&, const Domain&, int );


/******************************************************************
* A fast distance calculation                                     *
******************************************************************/
bool CalcDist3DIteration( Array<float> &Distance, const Domain & )
{
    const float sq2 = sqrt(2.0f);
    const float sq3 = sqrt(3.0f);
    float dist0[27];
    dist0[0] = sq3;   dist0[1] = sq2;   dist0[2] = sq3;
    dist0[3] = sq2;   dist0[4] = 1;     dist0[5] = sq2;
    dist0[6] = sq3;   dist0[7] = sq2;   dist0[8] = sq3;
    dist0[9] = sq2;   dist0[10] = 1;    dist0[11] = sq2;
    dist0[12] = 1;    dist0[13] = 0;    dist0[14] = 1;
    dist0[15] = sq2;  dist0[16] = 1;    dist0[17] = sq2;
    dist0[18] = sq3;  dist0[19] = sq2;  dist0[20] = sq3;
    dist0[21] = sq2;  dist0[22] = 1;    dist0[23] = sq2;
    dist0[24] = sq3;  dist0[25] = sq2;  dist0[26] = sq3;
    bool changed = false;
    for (size_t k=1; k<Distance.size(2)-1; k++) {
        for (size_t j=1; j<Distance.size(1)-1; j++) {
            for (size_t i=1; i<Distance.size(0)-1; i++) {
                float dist[27];
                dist[0] = Distance(i-1,j-1,k-1);   dist[1] = Distance(i,j-1,k-1);   dist[2] = Distance(i+1,j-1,k-1);
                dist[3] = Distance(i-1,j,k-1);     dist[4] = Distance(i,j,k-1);     dist[5] = Distance(i+1,j,k-1);
                dist[6] = Distance(i-1,j+1,k-1);   dist[7] = Distance(i,j+1,k-1);   dist[8] = Distance(i+1,j+1,k-1);
                dist[9] = Distance(i-1,j-1,k);     dist[10] = Distance(i,j-1,k);    dist[11] = Distance(i+1,j-1,k);
                dist[12] = Distance(i-1,j,k);      dist[13] = Distance(i,j,k);      dist[14] = Distance(i+1,j,k);
                dist[15] = Distance(i-1,j+1,k);    dist[16] = Distance(i,j+1,k);    dist[17] = Distance(i+1,j+1,k);
                dist[18] = Distance(i-1,j-1,k+1);  dist[19] = Distance(i,j-1,k+1);  dist[20] = Distance(i+1,j-1,k+1);
                dist[21] = Distance(i-1,j,k+1);    dist[22] = Distance(i,j,k+1);    dist[23] = Distance(i+1,j,k+1);
                dist[24] = Distance(i-1,j+1,k+1);  dist[25] = Distance(i,j+1,k+1);  dist[26] = Distance(i+1,j+1,k+1);
                float tmp = 1e100;
                for (int k=0; k<27; k++)
                    tmp = std::min(tmp,dist[k]+dist0[k]);
                if ( tmp < Distance(i,j,k) ) {
                    Distance(i,j,k) = tmp;
                    changed = true;
                }
            }
        }
    }
    return changed;
}
void CalcDist3D( Array<float> &Distance, const Array<char> &ID, const Domain &Dm )
{
    PROFILE_START("Calc Distance");
    // Initialize the distance to be 0 fore the cells adjacent to the interface
    Distance.fill(1e100);
    for (size_t k=1; k<ID.size(2)-1; k++) {
        for (size_t j=1; j<ID.size(1)-1; j++) {
            for (size_t i=1; i<ID.size(0)-1; i++) {
                char id = ID(i,j,k);
                if ( id!=ID(i-1,j,k) || id!=ID(i+1,j,k) || id!=ID(i,j-1,k) || id!=ID(i,j+1,k) || id!=ID(i,j,k-1) || id!=ID(i,j,k+1) )
                    Distance(i,j,k) = 0.5;
            }
        }
    }
    // Compute the distance everywhere
    fillHalo<float> fillData(Dm.Comm, Dm.rank_info,Dm.Nx,Dm.Ny,Dm.Nz,1,1,1,0,1);
    while ( true ) {
        // Communicate the halo of values
        fillData.fill(Distance);
        // The distance of the cell is the minimum of the distance of the neighbors plus the distance to that node
        bool changed = CalcDist3DIteration( Distance, Dm );
        changed = sumReduce(Dm.Comm,changed);
        if ( !changed )
            break;
    }
    // Update the sign of the distance
    for (size_t i=0; i<ID.length(); i++)
        Distance(i) *= ID(i)>0 ? 1:-1;
    PROFILE_STOP("Calc Distance");
}


/******************************************************************
* A fast distance calculation                                     *
******************************************************************/
void CalcDistMultiLevelHelper( Array<float> &Distance, const Domain &Dm )
{
    size_t ratio = 4;
    std::function<float(const Array<float>&)> coarsen = [ratio]( const Array<float>& data )
    {
        float tmp = 1e100;
        int nx = data.size(0);
        int ny = data.size(1);
        int nz = data.size(2);
        for (int k=0; k<nz; k++) {
            float z = k-0.5*(nz-1);
            for (int j=0; j<ny; j++) {
                float y = j-0.5*(ny-1);
                for (int i=0; i<nx; i++) {
                    float x = i-0.5*(nx-1);
                    tmp = std::min(data(i,j,k)+sqrt(x*x+y*y+z*z),tmp);
                }
            }
        }
        return tmp/ratio;
    };
    int Nx = Dm.Nx-2;
    int Ny = Dm.Ny-2;
    int Nz = Dm.Nz-2;
    ASSERT(int(Distance.size(0))==Nx+2&&int(Distance.size(1))==Ny+2&&int(Distance.size(2))==Nz+2);
    fillHalo<float> fillData(Dm.Comm,Dm.rank_info,Nx,Ny,Nz,1,1,1,0,1);
    if ( Nx%ratio==0 && Nx>8 && Ny%ratio==0 && Ny>8 && Nz%ratio==0 && Nz>8 ) {
        // Use recursive version
        int Nr = std::max(std::max(ratio,ratio),ratio);
        // Run Nr iterations, communicate, run Nr iterations
        for (int i=0; i<Nr; i++)
            CalcDist3DIteration( Distance, Dm );
        /*fillData.fill(Distance);
        for (int i=0; i<Nr; i++)
            CalcDist3DIteration( Distance, Dm );*/
        // Coarsen
        Array<float> dist(Nx,Ny,Nz);
        fillData.copy(Distance,dist);
        auto db = Dm.getDatabase()->cloneDatabase();
        auto n = db->getVector<int>( "n" );
        db->putVector<int>( "n", { n[0]/ratio, n[1]/ratio, n[2]/ratio } );
        Domain Dm2(db);
        Dm2.CommInit(Dm.Comm);
        fillHalo<float> fillData2(Dm2.Comm,Dm2.rank_info,Nx/ratio,Ny/ratio,Nz/ratio,1,1,1,0,1);
        auto dist2 = dist.coarsen( {ratio,ratio,ratio}, coarsen );
        Array<float> Distance2(Nx/ratio+2,Ny/ratio+2,Nz/ratio+2);
        fillData2.copy(dist2,Distance2);
        // Solve
        CalcDistMultiLevelHelper( Distance2, Dm2 );
        // Interpolate the coarse grid to the fine grid
        fillData2.copy(Distance2,dist2);
        for (int k=0; k<Nz; k++) {
            int k2 = k/ratio;
            float z = (k-k2*ratio)-0.5*(ratio-1);
            for (int j=0; j<Ny; j++) {
                int j2 = j/ratio;
                float y = (j-j2*ratio)-0.5*(ratio-1);
                for (int i=0; i<Nx; i++) {
                    int i2 = i/ratio;
                    float x = (i-i2*ratio)-0.5*(ratio-1);
                    dist(i,j,k) = std::min(dist(i,j,k),ratio*dist2(i2,j2,k2)+sqrt(x*x+y*y+z*z));
                }
            }
        }
        fillData.copy(dist,Distance);
        // Run Nr iterations, communicate, run Nr iterations
        for (int i=0; i<Nr; i++)
            CalcDist3DIteration( Distance, Dm );
        fillData.fill(Distance);
        for (int i=0; i<Nr; i++)
            CalcDist3DIteration( Distance, Dm );
    } else {
        // Use coarse-grid version
        while ( true ) {
            // Communicate the halo of values
            fillData.fill(Distance);
            // The distance of the cell is the minimum of the distance of the neighbors plus the distance to that node
            bool changed = CalcDist3DIteration( Distance, Dm );
            changed = sumReduce(Dm.Comm,changed);
            if ( !changed )
                break;
        }
    }
}
void CalcDistMultiLevel( Array<float> &Distance, const Array<char> &ID, const Domain &Dm )
{
    PROFILE_START("Calc Distance Multilevel");
    int Nx = Dm.Nx-2;
    int Ny = Dm.Ny-2;
    int Nz = Dm.Nz-2;
    ASSERT(int(Distance.size(0))==Nx+2&&int(Distance.size(1))==Ny+2&&int(Distance.size(2))==Nz+2);
    fillHalo<float> fillData(Dm.Comm,Dm.rank_info,Nx,Ny,Nz,1,1,1,0,1);
    // Initialize the distance to be 0 fore the cells adjacent to the interface
    Distance.fill(1e100);
    for (size_t k=1; k<ID.size(2)-1; k++) {
        for (size_t j=1; j<ID.size(1)-1; j++) {
            for (size_t i=1; i<ID.size(0)-1; i++) {
                char id = ID(i,j,k);
                if ( id!=ID(i-1,j,k) || id!=ID(i+1,j,k) || id!=ID(i,j-1,k) || id!=ID(i,j+1,k) || id!=ID(i,j,k-1) || id!=ID(i,j,k+1) )
                    Distance(i,j,k) = 0.5;
            }
        }
    }
    // Solve the for the distance using a recursive method
    CalcDistMultiLevelHelper( Distance, Dm );
    // Update the sign of the distance
    for (size_t i=0; i<ID.length(); i++)
        Distance(i) *= ID(i)>0 ? 1:-1;
    fillData.fill(Distance);
    // Run a quick filter to smooth the data
    float sigma = 0.6;
    Array<float> H = imfilter::create_filter<float>( { 1 }, "gaussian", &sigma );
    std::vector<imfilter::BC> BC(3,imfilter::BC::replicate);
    Distance = imfilter::imfilter_separable<float>( Distance, {H,H,H}, BC );
    PROFILE_STOP("Calc Distance Multilevel");
}
