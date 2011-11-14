function curves2chains_accordion

global thisCurve;

spinner = 1;
q00 = [0;0;0;0;0;0;0];
q00 = [-pi/3; -pi/3; -pi/3; 0; pi/3; pi/3; pi/3];

Nvp = 7;
N = 10;
Q = zeros(Nvp,N);

T = linspace(0,5,Nvp)';

Curves = struct('f',{@f1; @f1; @f4; @f2; @f5; @f3; @f3},...
               'dp',{@dp1; @dp1; @dp4; @dp2; @dp5; @dp3; @dp3},...
               'xw',{4.9; 1.9;   0.4;  2.6;  0.4;  2.; 4.9},...
               'q0',{[pi/3; zeros(N-1,1)];...
                     [pi/1.5; zeros(N-1,1)];...
                     cumsum([pi-1.9832    2.0501    3.2292    3.0540   -3.0540    3.0540   -3.0540    3.0540   -3.0540    2.0501])';...
                     cumsum([pi/2-1.9832    2.0501    3.2292    3.0540   -3.0540    3.0540   -3.0540    3.0540   -3.0540    2.0501])';...
                     cumsum([-1.9832    2.0501    3.2292    3.0540   -3.0540    3.0540   -3.0540    3.0540   -3.0540    2.0501])';...
                     [-pi/1.5; zeros(N-1,1)];...
                     [-pi/2.2; zeros(N-1,1)]});

close all;
h = figure;
hold on;
axis equal;     axis off;
for i = 1:Nvp
    % check length of function
    x = 0:.01: Curves(i).xw;
    plot( x, Curves(i).f(x) ,'r','LineWidth',i/7*2+1) % plot
    plot( Curves(i).xw , Curves(i).f( Curves(i).xw ) ,'r');
    text( Curves(i).xw , Curves(i).f( Curves(i).xw ) , sprintf(' curve %d',i) ,'VerticalAlignment','bottom');
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
        
%         Q(n,1)
%         Q(n,1) > pi
%         Q(n,1)
%         Q(n,1) = Q(n,1) - (Q(n,1) > pi)*2*pi
%         
%         if n == 3
%             Q(n,:) = [+pi/2-0.5230   -4.1448    3.1416    3.0524   -3.0524    3.0524   -3.0524    3.0524   -3.0524    2.0492];
%         elseif n == 5
%             Q(n,:) = [-pi/2+0.5230   -4.1448    3.1416    3.0524   -3.0524    3.0524   -3.0524    3.0524   -3.0524    2.0492];
%         end
    end

    Q2 = Q;
    
    for j = 2:N
        Q2(:,j) = Q(:,j) - Q(:,j-1);
    end

    Q2 = Q2 - ( Q2 >= pi )*2*pi + ( Q2 < -pi )*2*pi;
    
    figure(h);
    axis off;
    for i = 1:Nvp
        
        plotlink2( Q2(i,:) ); % plot links
        
    end
    
    if ~input('Did the curves succeed? Should I save?: ')
        return;
    else
        % convert to relative angles

        if spinner
            Q2 = [q00 Q2];
        end
        
%         Q2 = mod(Q2,2*pi);

        OUT = [Q2 T];
        save path4accordion.mat OUT -ascii
    end    
end

end % curves2chains

function f = f1(x)
f = 5*(1-exp(-5*x));
end
function f = f2(x)
f = 0*x;
end
function f = f4(x)
f = 9*x;
end
function f = f3(x)
f = -5*( 1-exp(-5*x) );
end
function f = f5(x)
f = -f4(x);
end
function dp = dp1(x)
df = 25*exp(-5*x);
dp = sqrt( 1 + df.^2 );
end
function dp = dp4(x)
[a b] = size(x);
df = 9*ones(a,b);
dp = sqrt( 1 + df.^2 );
end
function dp = dp2(x)
df = 0*x;
dp = sqrt( 1 + df.^2 );
end
function dp = dp3(x)
df = -25*exp(-5*x);
dp = sqrt( 1 + df.^2 );
end
function dp = dp5(x)
dp = dp4(x);
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


function plotlink2(q)

N = length(q);

x = zeros(1,N+1);
y = zeros(1,N+1);

for i = 2:N+1

    for j = 1:i-1
        x(i) = x(i) + cos(sum(q(1:j)));
        y(i) = y(i) + sin(sum(q(1:j)));
    end
    plot([x(i-1) x(i)],[y(i-1) y(i)])
    plot(x(i),y(i),'.')
    text(mean(x(i-1:i)),mean(y(i-1:i)),sprintf('%d',i-1))
end

end