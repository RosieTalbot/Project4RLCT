t =  -5:0.5:5;


for i = 1:length(t)



xt1 = @(t) t;
yt1 = @(t) cos(t);
zt1 = @(t) sin(t);
xt2 = @(t) t;
yt2 = @(t) cos(2*t);
zt2 = @(t) sin(2*t);
fplot3(xt1,yt1,zt1)
hold on
fplot3(xt2,yt2,zt2)

