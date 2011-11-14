function q = curve2chain

global xw;

xw = 12;

s = quadgk(@pathint,0,xw)

N = 15;

q0 = zeros(N,1);
q0(1) = pi/10;

[q,fval,exitflag] = fsolve(@fn,q0);

fprintf('exit flag is %d\n',exitflag);

plotlink(q);


end

function f = mycurve(x)

f = sin(x);
%f = x;

end

function f = pathint(x)

f = sin(x);

f = sqrt( 1 + (cos(x)).^2 );
%f = x;

end

function F = fn(Q)

global xw;

N = length(Q);

F = zeros(N,1);

yw = mycurve(xw);

x = zeros(N,1);
y = zeros(N,1);
x(1) = cos(Q(1));
y(1) = sin(Q(1));
for i = 2:length(Q)
    x(i) = x(i-1) + cos(Q(i));
	y(i) = y(i-1) + sin(Q(i));
end

for i = 1:N-2
    F(i) = mycurve( mean(x(i:i+1)) ) - mean(y(i:i+1));
end

F(N-1) = xw - x(N);
F(N) = yw - y(N);

end

function plotlink(q)

global xw;

N = length(q);

x = zeros(1,N+1);
y = zeros(1,N+1);

close all;
figure;
hold on;
axis equal;
for i = 2:N+1

    for j = 1:i-1
        x(i) = x(i) + cos(q(j));
        y(i) = y(i) + sin(q(j));
    end
    plot([x(i-1) x(i)],[y(i-1) y(i)])
    text(x(i),y(i),sprintf('%d',i-1))
end

xp = 0:.01:xw;
plot(xp,mycurve(xp),'r')

end