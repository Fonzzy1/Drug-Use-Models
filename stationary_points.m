syms S U A R p1 p2 p3 p4 p5 b ds du da dr




ode1 =  ((-S*p1*(U+A) + U*p2*(S+R) -ds*S + b*(S+U+A+R))*(S+U+A+R) - (b*(S + U + A + R) - ds*S -du*U - da*A - dr*R)*S)/(S+U+A+R)^2;
ode2 =  ((S*p1*(U+A) - U*p2*(S+R) -U*p3 -du*U)*(S+U+A+R) - (b*(S + U + A + R) - ds*S -du*U - da*A - dr*R)*U)/(S+U+A+R)^2;
ode3 =  ((U*p3 - A*p4*(S+R) + R*p5*(U+A) - da*A)*(S+U+A+R) - (b*(S + U + A + R) - ds*S -du*U - da*A - dr*R)*A)/(S+U+A+R)^2;
ode4 =  (( A*p4*(S+R) - R*p5*(U+A) - dr*R)*(S+U+A+R) - (b*(S + U + A + R) - ds*S -du*U - da*A - dr*R)*R)/(S+U+A+R)^2;

odes = [ode1, ode2, ode3, ode4];


j = jacobian(odes,[S,U,A,R]);


% Define New function for me to minimise
h = 0.1;
cols = ['S','U','A','R','ds','du','da','dr','mag'];

res = [];


for S = 0:h:1
    for U = 0:h:1
        for A = 0:h:1
            R = 1 - S - U - A;
            if R >= 0
                    
                    x = [S;U;A;R];
                    d = (subs(j) * x);
                    v = [S,U,A,R,d(1),d(2),d(3),d(4),sqrt(d(1)^2 + d(2)^2 +d(3)^2 +d(4)^2)] ;
                    res = [res;v];
            end
        end
    end
end

res

sympref('FloatingPointOutput',true)

[p1, p2, p3, p4, p5, b, ds, du, da, dr] = deal(0.01,0.001,0.001, 0.001, 0.001,0.001,0.0001,0.0005,0.001,0.0001);

subs(res)

1 + r+s_

