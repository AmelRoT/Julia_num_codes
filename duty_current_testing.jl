T =320
pwm_d = 250
duty = pwm_d /T;

I_in = 0.25
I_out = I_in*duty

I_measured = 0.168
delta_I = I_out-I_measured

delta_I/I_out*100
