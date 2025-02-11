function PlotCircle(ordparam,U)
th = 0:pi/50:2*pi;
S = sin(U);
C = cos(U);

x = real(ordparam);
y = imag(ordparam);


p = plot(C,S,'ob','MarkerSize',10)
p.MarkerFaceColor = [0 0 1];
hold on
plot(cos(th),sin(th),'-r','LineWidth', 3)
quiver(0,0,x,y,0,'.m')
plot(0,0,'.k','MarkerSize',30)
plot(x,y,'.m','MarkerSize',40)
axis square
xlim([-1,1])
ylim([-1,1])
title('Angles of Particles at tMax')
legend('Particles','Unit Circle','','Origin','OrderParameter')
legend('off')
fontsize(24,"points")
end