T_2 = [0.4	0.05
0.7	0.25
-0.1	0.25
-0.35	0.52
0.25	0.5
1.13	0.4
1.12	0.81
1.3	-0.13
-0.72	0.55
0	0.9
0	0
0.3	0.11
0.25	0.33
0.56	0.3
0.6	1.01
0.95	0.6
];


%% main

GPSstd = 5;
GPSerrE = 0;
GPSerrN = 0.0;
wtHdg = 0; wtHdg = wtHdg * 180/pi;

e_g = atan(std(T_2)*2/407);
e_g_y = e_g(1) * 1;
e_g_p = e_g(2) * 0;

GPSerrX = GPSerrN * sin(wtHdg) + GPSerrE * cos(wtHdg);
GPSerrY = GPSerrN * cos(wtHdg) + GPSerrN * sin(wtHdg);
lerr = 20 + GPSerrY;

P_err = lerr*tan(e_g_y)+GPSerrX

P_err_angle = atan(P_err/20)*180/pi

P_err_percent = P_err / 1.04 * 100