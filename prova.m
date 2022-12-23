clc; clear; close all;

t = 0:0.1:60;
x = 20*t; y = cos(t); z = sin(t);
view(3);
hold on; grid on;
color = jet(size(t,2)-1);
for k = 2:size(t,2)
    plot3(x(k-1:k),y(k-1:k),z(k-1:k),'color',color(k-1,:))
    axis([0 1800.7 -1.7 1.7 -1.3 1.3])
    
    pause(0.01)
end
colorbar