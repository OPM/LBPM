//--------------------------------------------------------------------------------------------------------
inline double ContactAngle(DTPoint3D &pt, DTVectorField3D &gradF, DTVectorField3D &gradS)
{
	double theta;
	int m = gradS.Grid().m();
	int n = gradS.Grid().n();
	int o = gradS.Grid().o();

	double dx = gradS.Grid().dx();
	double dy = gradS.Grid().dy();
	double dz = gradS.Grid().dz();

	DTPoint3D origin = gradS.Grid().Origin();

	// Get the gradient arrays
	DTFloatArray Fx = gradF.X();
	DTFloatArray Fy = gradF.Y();
	DTFloatArray Fz = gradF.Z();

	DTFloatArray Sx = gradS.X();
	DTFloatArray Sy = gradS.Y();
	DTFloatArray Sz = gradS.Z();

	int i,j,k;
	// Determine the cube containing this point
//	i = int(floor((pt.x-origin.x)/dx));	// real space
//	j = int(floor((pt.y-origin.y)/dy));
//	k = int(floor((pt.z-origin.z)/dz));
	i = int(floor(pt.x));
	j = int(floor(pt.y));
	k = int(floor(pt.z));

	double x,y,z;	// map to 0,1 cube
//	x = (pt.x-origin.x-i*dx)/dx;  // real space
//	y = (pt.y-origin.y-j*dy)/dy;
//	z = (pt.z-origin.z-k*dz)/dz;
	x = pt.x - double(i);
	y = pt.y - double(j);
	z = pt.z - double(k);

	// evaluate contact angle at corners of the cube
	DTMutableDoubleArray cube(2,2,2);
	// theta = acos ( -(gradF*gradS) / (|gradF| |gradS|) )
	cube(0,0,0) = -( Fx(i,j,k)*Sx(i,j,k)+Fy(i,j,k)*Sy(i,j,k)+Fz(i,j,k)*Sz(i,j,k) )
					   /( sqrt(pow(Fx(i,j,k),2)+pow(Fy(i,j,k),2)+pow(Fz(i,j,k),2))
						  *sqrt(pow(Sx(i,j,k),2)+pow(Sy(i,j,k),2)+pow(Sz(i,j,k),2)) );
	cube(1,0,0) = -( Fx(i+1,j,k)*Sx(i+1,j,k)+Fy(i+1,j,k)*Sy(i+1,j,k)+Fz(i+1,j,k)*Sz(i+1,j,k) )
					   /( sqrt(pow(Fx(i+1,j,k),2)+pow(Fy(i+1,j,k),2)+pow(Fz(i+1,j,k),2))
						  *sqrt(pow(Sx(i+1,j,k),2)+pow(Sy(i+1,j,k),2)+pow(Sz(i+1,j,k),2)) );
	cube(0,1,0) = -( Fx(i,j+1,k)*Sx(i,j+1,k)+Fy(i,j+1,k)*Sy(i,j+1,k)+Fz(i,j+1,k)*Sz(i,j+1,k) )
					   /( sqrt(pow(Fx(i,j+1,k),2)+pow(Fy(i,j+1,k),2)+pow(Fz(i,j+1,k),2))
						  *sqrt(pow(Sx(i,j+1,k),2)+pow(Sy(i,j+1,k),2)+pow(Sz(i,j+1,k),2)) );
	cube(0,0,1) = -( Fx(i,j,k+1)*Sx(i,j,k+1)+Fy(i,j,k+1)*Sy(i,j,k+1)+Fz(i,j,k+1)*Sz(i,j,k+1) )
					   /( sqrt(pow(Fx(i,j,k+1),2)+pow(Fy(i,j,k+1),2)+pow(Fz(i,j,k+1),2))
						  *sqrt(pow(Sx(i,j,k+1),2)+pow(Sy(i,j,k+1),2)+pow(Sz(i,j,k+1),2)) );
	cube(1,1,0) = -( Fx(i+1,j+1,k)*Sx(i+1,j+1,k)+Fy(i+1,j+1,k)*Sy(i+1,j+1,k)+Fz(i+1,j+1,k)*Sz(i+1,j+1,k) )
					   /( sqrt(pow(Fx(i+1,j+1,k),2)+pow(Fy(i+1,j+1,k),2)+pow(Fz(i+1,j+1,k),2))
						  *sqrt(pow(Sx(i+1,j+1,k),2)+pow(Sy(i+1,j+1,k),2)+pow(Sz(i+1,j+1,k),2)) );
	cube(1,0,1) = -( Fx(i+1,j,k+1)*Sx(i+1,j,k+1)+Fy(i+1,j,k+1)*Sy(i+1,j,k+1)+Fz(i+1,j,k+1)*Sz(i+1,j,k+1) )
					   /( sqrt(pow(Fx(i+1,j,k+1),2)+pow(Fy(i+1,j,k+1),2)+pow(Fz(i+1,j,k+1),2))
						  *sqrt(pow(Sx(i+1,j,k+1),2)+pow(Sy(i+1,j,k+1),2)+pow(Sz(i+1,j,k+1),2)) );
	cube(0,1,1) = -( Fx(i,j+1,k+1)*Sx(i,j+1,k+1)+Fy(i,j+1,k+1)*Sy(i,j+1,k+1)+Fz(i,j+1,k+1)*Sz(i,j+1,k+1) )
					   /( sqrt(pow(Fx(i,j+1,k+1),2)+pow(Fy(i,j+1,k+1),2)+pow(Fz(i,j+1,k+1),2))
						  *sqrt(pow(Sx(i,j+1,k+1),2)+pow(Sy(i,j+1,k+1),2)+pow(Sz(i,j+1,k+1),2)) );
	cube(1,1,1) = -( Fx(i+1,j+1,k+1)*Sx(i+1,j+1,k+1)+Fy(i+1,j+1,k+1)*Sy(i+1,j+1,k+1)+Fz(i+1,j+1,k+1)*Sz(i+1,j+1,k+1) )
					   /( sqrt(pow(Fx(i+1,j+1,k+1),2)+pow(Fy(i+1,j+1,k+1),2)+pow(Fz(i+1,j+1,k+1),2))
						  *sqrt(pow(Sx(i+1,j+1,k+1),2)+pow(Sy(i+1,j+1,k+1),2)+pow(Sz(i+1,j+1,k+1),2)) );

	// construct polynomial approximation for contact angle
	double a,b,c,d,e,f,g,h; // trilinear coefficients: f(x,y,z) = a+bx+cy+dz+exy+fxz+gyz+hxyz
	a = cube(0,0,0);
	b = cube(1,0,0)-a;
	c = cube(0,1,0)-a;
	d = cube(0,0,1)-a;
	e = cube(1,1,0)-a-b-c;
	f = cube(1,0,1)-a-b-d;
	g = cube(0,1,1)-a-c-d;
	h = cube(1,1,1)-a-b-c-d-e-f-g;

	// evaluate at x,y,z
	theta = acos (a+b*x+c*y+d*z+e*x*y+f*x*z+g*y*z+h*x*y*z );

	return theta;
}
