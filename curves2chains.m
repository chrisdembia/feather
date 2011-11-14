function curves2chains

global thisCurve;

spinner = 1;
q00 = [0;0;0];

Nvp = 3;
N = 10;
Q = zeros(Nvp,N);

T = [0; 1; 2];

Curves = struct('f',{@f1; @f2; @f3},...
               'dp',{@dp1; @dp2; @dp3},...
               'xw',{4.9; 2.2; 4.9},...
               'q0',{[pi/2; zeros(N-1,1)];...
                      % 1      2       3  4      5       6      7    8  9 10      
                     %[pi/2.5 pi/10 -pi/5 pi/5 -pi/1.1 pi/2.5 -pi/1.5 pi/5 -pi/2.5 pi/2.5];...
                     [pi/2.1 0 0 0 0 -pi/1.2 0 0 0 0];...
                     [-pi/2; zeros(N-1,1)]});

close all;
h = figure;
hold on;
axis equal;
for i = 1:Nvp
    % check length of function
    x = 0:.01: Curves(i).xw;
    plot( x, Curves(i).f(x) ,'r') % plot
    plot( Curves(i).xw , Curves(i).f( Curves(i).xw ) ,'r');
    text( Curves(i).xw , Curves(i).f( Curves(i).xw ) , sprintf('curve %d',i) );
    s = quadgk( Curves(i).dp , 0, Curves(i).xw );
    fprintf('Path length of curve %d: %f\n', i, s);
end

if ~input('Continue? [1 for yes, 0 for no]: ')
    return;
else
    
    for n = 1:Nvp
        
        thisCurve = Curves(n);
        
        if n == 2
            [q,fval,exitflag] = fsolve( @chaineqns , Q(1,:)' );
        else
            [q,fval,exitflag] = fsolve( @chaineqns , thisCurve.q0 );
        end
        fprintf('Exit flag for curve %d is %d\n',n,exitflag);
        
        Q(n,:) = q';
    end

    figure(h);
    for i = 1:Nvp
        plotlink( Q(i,:) ); % plot links
    end
    
    if ~input('Did the curves succeed? Should I save?: ')
        return;
    else
        % convert to relative angles
        Q2 = Q;
        for j = 2:N
            Q2(:,j) = Q(:,j) - Q(:,j-1);
        end
        
        if spinner
            Q2 = [q00 Q2];
        end
        
        OUT = [Q2 T];
        save path1.mat OUT -ascii
    end    
end

end % curves2chains

function f = f1(x)

%f = 0.1*x.^2;
%f = 0*x;
f = 5*(1-exp(-5*x));

end

function f = f2(x)

% f = 2*sqrt(x);
xo = 2.2;
ym = 4.5;

c = [xo^2 xo; xo^2/4 xo/2]\[0; ym];

f = c(1) * x.^2 + c(2) * x;

end

function f = f3(x)

f = -5*( 1-exp(-5*x) );

end

function dp = dp1(x)

% df = 2*.1*x;
% df = 0*x;
df = 25*exp(-5*x);
dp = sqrt( 1 + df.^2 );

end

function dp = dp2(x)

xo = 2.2;
ym = 4.5;

c = [xo^2 xo; xo^2/4 xo/2]\[0; ym];

df = 2* c(1) * x + c(2);
% df = 2*0.5./sqrt(x);
dp = sqrt( 1 + df.^2 );

end

function dp = dp3(x)

% df = -8*.2*exp(-.2*x);
df = -25*exp(-5*x);
dp = sqrt( 1 + df.^2 );

end

function F = chaineqns(Q)

global thisCurve N;

xw = thisCurve.xw;

N = length(Q);

F = zeros(N,1);

yw = thisCurve.f(xw);

x = zeros(N,1);
y = zeros(N,1);
x(1) = cos(Q(1));
y(1) = sin(Q(1));
for i = 2:length(Q)
    x(i) = x(i-1) + cos(Q(i));
	y(i) = y(i-1) + sin(Q(i));
end

for i = 1:N-2
    F(i) = thisCurve.f( mean(x(i:i+1)) ) - mean(y(i:i+1));
end

F(N-1) = xw - x(N);
F(N) = yw - y(N);

end

function plotlink(q)

global xw;

N = length(q);

x = zeros(1,N+1);
y = zeros(1,N+1);

for i = 2:N+1

    for j = 1:i-1
        x(i) = x(i) + cos(q(j));
        y(i) = y(i) + sin(q(j));
    end
    plot([x(i-1) x(i)],[y(i-1) y(i)])
    plot(x(i),y(i),'.')
    text(mean(x(i-1:i)),mean(y(i-1:i)),sprintf('%d',i-1))
end

end