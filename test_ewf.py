w_guess = 50000
WP = 10800
w_fuel = 0.3773

def calculate_ewf(w_guess):
    ewf = 0.93 * (w_guess ** -0.07)
    return ewf


while True:
    ewf = calculate_ewf(w_guess=w_guess)
    WTO_calc = WP / (1 - w_fuel - ewf)

    if abs(WTO_calc - w_guess) < 0.001:
        break

    w_guess = WTO_calc
    print(w_guess)

print(w_guess)