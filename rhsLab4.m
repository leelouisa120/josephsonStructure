function yp = rhsLab4(theta,t)

global I_o w gamma;

thetadot=I_o-cos(theta) +gamma*sin(w*t);

yp=thetadot;
end