




simpsons_weights = Dict(
    (:classic, 3)   => (1/3) .* (1, 4, 1),                                            # Simpson's Rule
    (:classic, 4)   => (3/8) .* (1, 3, 3, 1),                                         # Simpson's 3/8 Rule
    (:lownoise, 5) => (4/105) .* (11, 26, 31, 26, 11),                               # Low Noise Simpson's Rule Coefficients I
    (:lownoise, 6) => (5/336) .* (31, 61, 76, 76, 61, 31),                           # Low Noise Simpson's Rule Coefficients II
    (:lownoise, 7) => (1/14) .* (7, 12, 15, 16, 15, 12, 7)                           # Low Noise Simpson's Rule Coefficients III
)

booles_weights = Dict(
    (:classic, 5) => (2/45) .* (7, 32, 12, 32, 7),                                    # Boole's Rule I
    (:classic, 6) => (5/288) .* (19, 75, 50, 50, 75, 19),                             # Boole's Rule II 
    (:lownoise, 7) => (1/770) .* (268, 933, 786, 646, 786, 933, 268),                # Low Noise Boole's Rule Coefficients I
    (:lownoise, 8) => (7/31680) .* (1657, 5157, 4947, 4079, 4079, 4947, 5157, 1657), # Low Noise Boole's Rule Coefficients II
    (:lownoise, 9) => (8/6435) .* (309, 869, 904, 779, 713, 779, 904, 869, 309)      # Low Noise Boole's Rule Coefficients III
)