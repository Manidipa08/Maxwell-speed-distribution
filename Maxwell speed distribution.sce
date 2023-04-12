clc
clear
clf
//speed distribution
N = 2001
kb = 1.38D-23
R = 8.324
Na = 6.022D+23
h = 6.626D-34
v = linspace (0, 2000, N)
T = 601
m = 40*1.67*(1D-27)

function y=f(v)
    y = ((4*%pi)*(m/(2*%pi*kb*T))^(3/2)*v^2).*exp((-m*v^2)/(2*kb*T))
endfunction
d = f(v)
b(1:N) = v(1:N).*f(v(1:N))
c(1:N) = v(1:N)^2.*f(v(1:N))
avg = inttrap(v,b)
rms = sqrt(inttrap(v,c))
[ind , mp] = max(d)
theo_avg = sqrt((8*kb*T)/(%pi*m))
theo_rms = sqrt((3*kb*T)/m)
theo_mp = sqrt((2*kb*T)/(m))
show_window(1)
plot(v,f(v),'k')
plot(avg*ones(1,N),d,'-*-g')
plot(rms*ones(1,N),d,'-')
plot(mp*ones(1,N),d,'--c')
title('Maxwell Distribution Function (T = 600 K)')
legend('Function (v)','Average Velocity','RMS Value ','Most Probable Value')
xlabel("Velocity (m/s)")
ylabel('Function (v)')
xgrid()

disp("Theoretical value of Average Velocity: "+string(theo_avg))
disp("Calculated value of Average Velocity: "+string(avg))
disp("Theoretical value of RMS Velocity: "+string(theo_rms))
disp("Calculated value of RMS Velocity: "+string(rms))
disp("Theoretical value of Most Probable Velocity: "+string(theo_mp))
disp("Calculated value of Most Probable Velocity: "+string(v(mp)))

//different temperature-------------
temp = [300, 500, 1200, 1500, 2000, 2500]
for i = 1:6
    function y=f(v)
        y = ((4*%pi)*(m/(2*%pi*kb*temp(i)))^(3/2)*v^2).*exp((-m*v^2)/(2*kb*temp(i)))
    endfunction
    func = f(v)
    show_window(2)
    plot(v,func)
    title('Maxwell Distribution Function','Fontsize','5')
    legend('300 K','500 K','1200 K','1500 K','2000 K', '2500K')
    xlabel("Velocity (m/s)")
    ylabel('Function (v)')
end
//different gases--------------
m_atomic = [4 20 40 83 132 222]*1.67*(1D-27)
velo = linspace (0, 5000, N)
for j = 1:6
    function y=f(velo) 
        y = ((4*%pi)*(m_atomic(j)/(2*%pi*kb*T))^(3/2)*velo^2).*exp((-m_atomic(j)*velo^2)/(2*kb*T))
    endfunction
    func_m = f(velo)
    show_window(3)
    plot(velo,func_m)
    title('Maxwell Distribution Function (T = 600 K)')
    legend('He (4)','Ne (20)','Ar (40)','Kr (83)','Xe (132)', 'Rn (222)')
    xlabel("Velocity (m/s)")
    ylabel('Function (v)')
    xgrid()
end
//velocity distribution----------------
v_d = linspace (-1500, 1500, N)
function y=f(v_d) 
    y = ((m/(2*%pi*kb*T))^(1/2).*exp((-m*v_d^2)/(2*kb*T))) 
endfunction
q = f(v_d)
n(1:N) = v_d(1:N).*f(v_d(1:N))
l(1:N) = v_d(1:N)^2.*f(v_d(1:N))
av_d = inttrap(v_d,n)
rms_d = sqrt(inttrap(v_d,l))
[xx,mp_d] = max(q)
th_avd = 0 
th_rmsd = sqrt((kb*T)/m)
th_mpd = 0 
show_window(4)
plot(v_d,q,'r')
plot(rms_d*ones(1,N),q,'-*g','thickness',1.5)
plot(v_d(mp_d)*ones(1,N),q,'-o-c','thickness',2)
plot(av_d*ones(1,N),q,'--k','thickness',1.5)
title('Maxwell Distribution Function (T = 600 K)')
legend('Function (v)','Average Velocity','RMS','Most Probable')
xlabel("Velocity (m/s)")
ylabel('Function (v)')
xgrid()
disp("Theoretical value of Average Velocity: "+string(th_avd))
disp("Calculated value of Average Velocity: "+string(av_d))
disp("Theoretical value of RMS Velocity: "+string(th_rmsd))
disp("Calculated value of RMS Velocity: "+string(rms_d))
disp("Theoretical value of Most Probable Velocity: "+string(th_mpd))
disp("Calculated value of Most Probable Velocity: "+string(v_d(mp_d)))
//different temperature in velocity distribution-------------
v_dd = linspace (-2500, 2500, N)
for ii = 1:6
    function y=f(v_dd) 
        y = ((m/(2*%pi*kb*temp(ii)))^(1/2).*exp((-m*v_dd^2)/(2*kb*temp(ii)))), 
    endfunction
    func_T = f(1:N)
    show_window(5)
    plot(v_dd,f(v_dd))
    title('Maxwell Distribution Function ')
    legend('300 K','500 K','1200 K','1500 K','2000 K', '2500K')
    xlabel("Velocity (m/s)")
    ylabel('Function (v_dist)')
    xgrid()
end
//different gases for velocity distribution------------------
for jj = 1:6
    function y=f(v_dd)
        y = ((m_atomic(jj)/(2*%pi*kb*T))^(1/2).*exp((-m_atomic(jj)*v_dd^2)/(2*kb*T)))
     endfunction
    func_g = f(v_dd)
    show_window(6)
    plot(v_dd,func_g)
    title('Maxwell Distribution Function (T = 600 K)')
    legend('He (4)','Ne (20)','Ar (40)','Kr (83)','Xe (132)', 'Rn (222)')
    xlabel("Velocity (m/s)")
    ylabel('Function (v)')
    xgrid()
end
