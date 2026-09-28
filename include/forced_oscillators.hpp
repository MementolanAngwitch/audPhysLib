

//Forced oscillations, (w = omega, j is the complex constant)
struct Drive { double F; double w; };
//complex displacement x = (1/jw) * Fe^{jwt}/(Rm + j(wm-s/w))

//complex speed u = (Fe^{jwt})/(Rm + j(wm - s/w))

//impedance Zm = Rm + jXm

//reactance Xm = wm - s/w

//impedance magnitude Zm = (Rm^2 + (wm -s/w)^2)^(1/2)

//phase angle theta = tan^-1 (Xm/Rm) = tan-1 ((wm-s/w)/Rm) use atan2

//Zm = Zm e^(j theta)

//Simplfiied impedance Zm = f/u, f = Fexp(jwt)

//Simplified equations of motion
// u = f/Zm
// x = f/(j*w*Zm)

//Actual (real) value
// x = (F/(w*Zm) sin(wt - theta))
// u = (F/Zm)cos(wt - theta)

// Power Relations
//instantanetous power (watts) PI_i  = (F^2/Zm) cos(wt) * cos(wt - theta)
//average power (watts)        PI = F^2/(2*Zm) * cos(theta) = (F^2Rm)/(2Zm^2) 

//Mechanical Resonance 
//w0 is where Xm vanishes and Zm = Rm
//Resonance values
// u_res = (F/Rm) cos(w0 t)
// x_res = (F/(w0*Rm)) sin(w0 * t)
// Quality factor Q = w0/(w_u - w_l)
// w_u and w_l are the two angular frequenices above and below resoance at which avg is half its resonance value
// Q = w0m /Rm
// Q = w0/(2beta)
// Q = 1/2 w0 t (t = relaxation time)

//Linear Combination
//Same angular frequencies
// x = A cos(wt+phi)
// A = [(SUM(An*cos(phi_n))^2)^2 + (SUM(A_n*sin(phi_n)))^2]^(1/2)
// tan(phi) = SUM(A_n * sin(phi_n)) / SUM(A_n * cos(phi_n))
// different angular frequenices
// w2 = w1 + delta(w)
// x = Ae^(j(w1*t + phi))
// A = [A1^2 +A2^2 + 2A1A2 * cos(phi1 - phi2 - delta*w*t)]^1/2
// tan(phi) = (A1 * sin(phi_1) + A2 * sin(phi_2 + delta*w*t))/(A1 * cos(phi_1) + A2 cos(phi_2 + delta*w*t))
} 