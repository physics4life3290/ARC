# ================================================
# Fixed Implicit Euler solver for full PP chain
# ================================================

using LinearAlgebra
using Printf

# Net stellar energy deposition (MeV, neutrinos removed)

const Q_pp      = 0.420      # p+p (most energy lost to ν)
const Q_dp      = 5.494
const Q_33      = 12.860
const Q_3He7    = 1.586
const Q_7Be_e   = 1.0        # ~90% lost to ν
const Q_7Be8    = 0.137
const Q_8_decay = 6.6        # neutrino-corrected

# CNO (neutrino-corrected)
const Q_c12p    = 1.944
const Q_n13dec  = 1.2
const Q_c13p    = 7.551
const Q_n14p    = 7.297
const Q_o15dec  = 1.0
const Q_n15p    = 4.966


# -------------------------------
# Constants & Helpers
# -------------------------------
const mu = 1.660539e-24 

function reaclib_rate(T9, a)
    # 7-parameter REACLIB fit
    return exp(a[1] + a[2]/T9 + a[3]/T9^(1/3) + a[4]*T9^(1/3) +
               a[5]*T9 + a[6]*T9^(5/3) + a[7]*log(T9))
end

# -------------------------------
# Species Indexing
# -------------------------------
# 1:p, 2:d, 3:He3, 4:He4, 5:Be7, 6:B8
# 7:C12, 8:N13, 9:C13, 10:N14, 11:O15, 12:N15
# -------------------------------

function cno_pp_f(Y, rho, T9, coeffs)
    f = zeros(12)
    # Unpack coefficients (Simplified for brevity)
    (a_pp, a_dp, a_33, a_3He7, a_7Be8, a_7Be_e, a_8,
     a_c12p, a_n13_dec, a_c13p, a_n14p, a_o15_dec, a_n15p) = coeffs

    # Rates
    λ_pp     = 0.5 * reaclib_rate(T9, a_pp) * rho / mu
    λ_dp     = reaclib_rate(T9, a_dp) * rho / mu
    λ_33     = 0.5 * reaclib_rate(T9, a_33) * rho / mu
    λ_3He7   = reaclib_rate(T9, a_3He7) * rho / mu
    λ_c12p   = reaclib_rate(T9, a_c12p) * rho / mu
    λ_n13dec = reaclib_rate(T9, a_n13_dec)
    λ_c13p   = reaclib_rate(T9, a_c13p) * rho / mu
    λ_n14p   = reaclib_rate(T9, a_n14p) * rho / mu
    λ_o15dec = reaclib_rate(T9, a_o15_dec)
    λ_n15p   = reaclib_rate(T9, a_n15p) * rho / mu

    # Right Hand Side Equations
    # Proton consumption (PP + CNO)
    f[1] = -2*λ_pp*Y[1]^2 - λ_dp*Y[1]*Y[2] + 2*λ_33*Y[3]^2 - λ_c12p*Y[7]*Y[1] - λ_c13p*Y[9]*Y[1] - λ_n14p*Y[10]*Y[1] - λ_n15p*Y[12]*Y[1]
    
    # PP Chain results (He4 produced by n15p reaction too)
    f[4] = λ_33*Y[3]^2 + λ_n15p*Y[12]*Y[1] 

    # CNO Cycle species
    f[7]  = -λ_c12p*Y[7]*Y[1]  + λ_n15p*Y[12]*Y[1]   # C12: Consumed by p, produced by N15+p
    f[8]  =  λ_c12p*Y[7]*Y[1]  - λ_n13dec*Y[8]       # N13: Produced by C12+p, decays
    f[9]  =  λ_n13dec*Y[8]     - λ_c13p*Y[9]*Y[1]    # C13: Produced by N13 decay, consumed by p
    f[10] =  λ_c13p*Y[9]*Y[1]  - λ_n14p*Y[10]*Y[1]   # N14: Produced by C13+p, consumed by p
    f[11] =  λ_n14p*Y[10]*Y[1] - λ_o15dec*Y[11]      # O15: Produced by N14+p, decays
    f[12] =  λ_o15dec*Y[11]    - λ_n15p*Y[12]*Y[1]   # N15: Produced by O15 decay, consumed by p

    return f
end

# -------------------------------
# Analytical Jacobian
# -------------------------------
function cno_pp_jacobian(Y, rho, T9, coeffs)
    # Unpack Abundances
    # 1:p, 2:d, 3:He3, 4:He4, 5:Be7, 6:B8
    # 7:C12, 8:N13, 9:C13, 10:N14, 11:O15, 12:N15
    Yp, Yd, Y3, Y4, Y7, Y8, Yc12, Yn13, Yc13, Yn14, Yo15, Yn15 = Y

    # Unpack Coeffs
    (a_pp, a_dp, a_33, a_3He7, a_7Be8, a_7Be_e, a_8,
     a_c12p, a_n13dec, a_c13p, a_n14p, a_o15dec, a_n15p) = coeffs

    # Rate Constants
    k_pp      = 0.5 * reaclib_rate(T9, a_pp) * rho / mu
    k_dp      =       reaclib_rate(T9, a_dp) * rho / mu
    k_33      = 0.5 * reaclib_rate(T9, a_33) * rho / mu
    k_3He7    =       reaclib_rate(T9, a_3He7) * rho / mu
    k_7Be8    =       reaclib_rate(T9, a_7Be8) * rho / mu
    k_7Be_e   =       reaclib_rate(T9, a_7Be_e)
    k_8       =       reaclib_rate(T9, a_8)
    
    k_c12p    =       reaclib_rate(T9, a_c12p) * rho / mu
    k_n13dec  =       reaclib_rate(T9, a_n13dec)
    k_c13p    =       reaclib_rate(T9, a_c13p) * rho / mu
    k_n14p    =       reaclib_rate(T9, a_n14p) * rho / mu
    k_o15dec  =       reaclib_rate(T9, a_o15dec)
    k_n15p    =       reaclib_rate(T9, a_n15p) * rho / mu

    J = zeros(12, 12)

    # --- 1. Protons (Y1) ---
    J[1,1]  = -4*k_pp*Yp - k_dp*Yd - k_7Be8*Y7 - k_c12p*Yc12 - k_c13p*Yc13 - k_n14p*Yn14 - k_n15p*Yn15
    J[1,2]  = -k_dp*Yp
    J[1,3]  = 4*k_33*Y3
    J[1,5]  = -k_7Be8*Yp
    J[1,7]  = -k_c12p*Yp
    J[1,9]  = -k_c13p*Yp
    J[1,10] = -k_n14p*Yp
    J[1,12] = -k_n15p*Yp

    # --- 2. Deuterium (Y2) ---
    J[2,1] = 2*k_pp*Yp - k_dp*Yd
    J[2,2] = -k_dp*Yp

    # --- 3. He3 (Y3) ---
    J[3,1] = k_dp*Yd
    J[3,2] = k_dp*Yp
    J[3,3] = -4*k_33*Y3 - k_3He7*Y4
    J[3,4] = -k_3He7*Y3

    # --- 4. He4 (Y4) ---
    J[4,1] = k_n15p*Yn15
    J[4,3] = 2*k_33*Y3
    J[4,5] = 2*k_7Be_e
    J[4,6] = 2*k_8
    J[4,12]= k_n15p*Yp

    # --- 5. Be7 (Y5) ---
    J[5,1] = -k_7Be8*Y7
    J[5,3] = k_3He7*Y4
    J[5,4] = k_3He7*Y3
    J[5,5] = -k_7Be8*Yp - k_7Be_e

    # --- 6. B8 (Y6) ---
    J[6,1] = k_7Be8*Y7
    J[6,5] = k_7Be8*Yp
    J[6,6] = -k_8

    # --- 7. C12 (Y7) ---
    J[7,1] = -k_c12p*Yc12 + k_n15p*Yn15
    J[7,7] = -k_c12p*Yp
    J[7,12]=  k_n15p*Yp

    # --- 8. N13 (Y8) ---
    J[8,1] = k_c12p*Yc12
    J[8,7] = k_c12p*Yp
    J[8,8] = -k_n13dec

    # --- 9. C13 (Y9) ---
    J[9,1] = -k_c13p*Yc13
    J[9,8] = k_n13dec
    J[9,9] = -k_c13p*Yp

    # --- 10. N14 (Y10) ---
    J[10,1] = k_c13p*Yc13 - k_n14p*Yn14
    J[10,9] = k_c13p*Yp
    J[10,10]= -k_n14p*Yp

    # --- 11. O15 (Y11) ---
    J[11,1] = k_n14p*Yn14
    J[11,10]= k_n14p*Yp
    J[11,11]= -k_o15dec

    # --- 12. N15 (Y12) ---
    J[12,1] = -k_n15p*Yn15
    J[12,11]= k_o15dec
    J[12,12]= -k_n15p*Yp

    return J
end

# Mass excesses in MeV (simplified values)
# Formula: Q = sum( dY/dt * (Mass_Excess_i) ) * conversion_factors
const mass_excess = [
    7.2890,  # p
    13.1357, # d
    14.9312, # He3
    2.4249,  # He4
    15.769,  # Be7
    22.921,  # B8
    0.0,     # C12 (Reference)
    5.345,   # N13
    3.125,   # C13
    2.863,   # N14
    2.855,   # O15
    0.101    # N15
]

function energy_generation(Y, rho, T9, coeffs)
    MeV_to_erg = 1.60218e-6
    NA = 6.02214076e23

    (a_pp, a_dp, a_33, a_3He7, a_7Be8, a_7Be_e, a_8,
     a_c12p, a_n13dec, a_c13p, a_n14p, a_o15dec, a_n15p) = coeffs

    # Reaction rates per gram per second
    r_pp    = 0.5 * reaclib_rate(T9, a_pp) * rho / mu * Y[1]^2
    r_dp    =       reaclib_rate(T9, a_dp) * rho / mu * Y[1]*Y[2]
    r_33    = 0.5 * reaclib_rate(T9, a_33) * rho / mu * Y[3]^2
    r_3He7  =       reaclib_rate(T9, a_3He7) * rho / mu * Y[3]*Y[4]
    r_7Be8  =       reaclib_rate(T9, a_7Be8) * rho / mu * Y[5]*Y[1]
    r_7Be_e =       reaclib_rate(T9, a_7Be_e) * Y[5]
    r_8     =       reaclib_rate(T9, a_8) * Y[6]

    r_c12p  =       reaclib_rate(T9, a_c12p) * rho / mu * Y[7]*Y[1]
    r_n13   =       reaclib_rate(T9, a_n13dec) * Y[8]
    r_c13p  =       reaclib_rate(T9, a_c13p) * rho / mu * Y[9]*Y[1]
    r_n14p  =       reaclib_rate(T9, a_n14p) * rho / mu * Y[10]*Y[1]
    r_o15   =       reaclib_rate(T9, a_o15dec) * Y[11]
    r_n15p  =       reaclib_rate(T9, a_n15p) * rho / mu * Y[12]*Y[1]

    Qdot = (
        Q_pp*r_pp +
        Q_dp*r_dp +
        Q_33*r_33 +
        Q_3He7*r_3He7 +
        Q_7Be8*r_7Be8 +
        Q_7Be_e*r_7Be_e +
        Q_8_decay*r_8 +
        Q_c12p*r_c12p +
        Q_n13dec*r_n13 +
        Q_c13p*r_c13p +
        Q_n14p*r_n14p +
        Q_o15dec*r_o15 +
        Q_n15p*r_n15p
    )

    return Qdot * NA * MeV_to_erg
end



# -------------------------------
# Updated Implicit Euler step for 12 species
# -------------------------------
# Atomic Mass (A) and Atomic Number (Z) for our 12 species
const A = [1, 2, 3, 4, 7, 8, 12, 13, 13, 14, 15, 15] # Mass numbers
const Z = [1, 1, 2, 2, 4, 5,  6,  7,  6,  7,  8,  7] # Proton numbers (Charge)

function implicit_euler_step!(Y, dt, rho, T9, coeffs; tol=1e-12, max_iter=100)
    Ynew = copy(Y)
    n = length(Y)
    In = Matrix{Float64}(I, n, n)
    
    # Calculate the target Baryon count and Charge count from initial state
    target_baryon = sum(Y .* A)
    target_charge = sum(Y .* Z)

    for k in 1:max_iter
        f = cno_pp_f(Ynew, rho, T9, coeffs)
        F = Ynew - Y - dt * f

        Jf = cno_pp_jacobian(Ynew, rho, T9, coeffs)
        J = In - dt * Jf

        ΔY = J \ F
        Ynew -= ΔY

        # --- CONSERVATION ENFORCEMENT (PROJECTION) ---
        # 1. Enforce Baryon Conservation: Normalize Mass Fraction
        mass_fraction_sum = sum(Ynew .* A)
        Ynew .*= (target_baryon / mass_fraction_sum)

        # 2. Convergence check
        if maximum(abs.(ΔY)) < tol
            return clamp.(Ynew, 0.0, Inf)
        end

        
    end

    eps = energy_generation(Ynew, rho, T9, coeffs)
    return clamp.(Ynew, 0.0, Inf), eps

end

# -------------------------------
# Execution & Initial Conditions
# -------------------------------

# 1:p, 2:d, 3:He3, 4:He4, 5:Be7, 6:B8
# 7:C12, 8:N13, 9:C13, 10:N14, 11:O15, 12:N15
Y = zeros(12)
# 1. Define desired Mass Fractions (X) - Must sum to 1.0
X_p   = 0.730   # 73% Hydrogen
X_he4 = 0.268   # 25% Helium
X_c12 = 0.001   # 1% Carbon catalyst
X_n14 = 0.001   # 1% Nitrogen catalyst

# 2. Convert Mass Fraction (X) to Molar Abundance (Y)
# Y_i = X_i / A_i
Y = zeros(12)
Y[1] = X_p / 1       # p
Y[4] = X_he4 / 4     # He4
Y[7] = X_c12 / 12    # C12
Y[10] = X_n14 / 14   # N14

# 3. Verify Initial Mass (Should be 1.0)
initial_mass = sum(Y .* A) 
println("Initial Mass Fraction Sum: ", initial_mass)

rho = 150.0          # Solar core density
T9 = 0.0175           # 20 Million K (CNO starts becoming visible here)

# --- Expanded Coeffs ---
# (Note: Using simplified placeholder exponents for the CNO rates)
a_pp      = [-43.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0]
a_dp      = [-15.0, -0.2, 0.0, 0.0, 0.0, 0.0, 0.0]
a_33      = [-15.0, -0.8, 0.0, 0.0, 0.0, 0.0, 0.0]
a_3He7    = [-18.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0]
a_7Be8    = [-20.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0]
a_7Be_e   = [-12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
a_8_decay = [log(2)/0.77, 0, 0, 0, 0, 0, 0]

# CNO Reactions
a_c12p    = [-25.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
a_n13dec  = [log(2)/597.0, 0, 0, 0, 0, 0, 0] # Half life ~10 mins
a_c13p    = [-20.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
a_n14p    = [-30.0, -1.2, 0.0, 0.0, 0.0, 0.0, 0.0] # Bottleneck
a_o15dec  = [log(2)/122.0, 0, 0, 0, 0, 0, 0] # Half life ~2 mins
a_n15p    = [-15.0, -0.8, 0.0, 0.0, 0.0, 0.0, 0.0]

coeffs = (
    a_pp, a_dp, a_33, a_3He7, a_7Be8, a_7Be_e, a_8_decay,
    a_c12p, a_n13dec, a_c13p, a_n14p, a_o15dec, a_n15p
)

# -------------------------------
# Simulation
# -------------------------------
dt = 1e13 # Large step to see the catalytic shift
Y_final, eps = implicit_euler_step!(Y, dt, rho, T9, coeffs)


# -------------------------------
# Output Results
# -------------------------------
names = ["p", "d", "He3", "He4", "Be7", "B8", "C12", "N13", "C13", "N14", "O15", "N15"]
# Atomic mass numbers for the 12 species


# Mass fraction sum (should be ~1.0)
sum_initial = sum(Y .* A)
sum_final = sum(Y_final .* A)

println("Initial Mass Fraction Sum: ", sum_initial)
println("Final Mass Fraction Sum:   ", sum_final)
@printf("\nEnergy generation rate: %.3e erg g⁻¹ s⁻¹\n", eps)


for i in 1:12
    @printf("%-5s: %.4e\n", names[i], Y_final[i])
end
