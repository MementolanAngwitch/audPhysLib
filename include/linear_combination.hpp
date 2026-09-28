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