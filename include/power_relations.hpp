#include "force_oscillators.hpp"

namespace audphys{
// Power Relations
//instantanetous power (watts) PI_i  = (F^2/Zm) cos(wt) * cos(wt - theta)
//average power (watts)        PI = F^2/(2*Zm) * cos(theta) = (F^2Rm)/(2Zm^2) 

inline double instant_power(const DampedOscillator& o, const Drive& d) {
	
}
}