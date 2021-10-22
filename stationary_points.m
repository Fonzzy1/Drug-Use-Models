syms s u a r p1 p2 p3 p4 p5 ds du da dr


ode1 =  (-s*p1*(u+a)+u*p2*(s+r)  -ds*s + 1)-s *(1 - ds*s - du*u - da*a - dr*r);
ode2 =  (s*p1*(u+a)-u*p2*(s+r) - u*p3 -du*u )-u*(1 - ds*s - du*u - da*a - dr*r) ;
ode3 =  (u*p3 -a*p4*(s+r) + r*p5*(u+a) -da*a )-a*(1 - ds*s - du*u - da*a - dr*r);
ode4 =  (a*p4*(s+r) - r*p5*(u+a) -dr*r )-r*(1 - ds*s -du*u - da*a - dr*r);
sums = s+u+a+r;



odes = [ode1 == 0, ode2 == 0, ode3 == 0 , ode4 == 0, sums  == 1, s>=0, u >= 0, a >=0, r>=0];


ss = [];
avg = [];

k = 3;

for p1i = -k:k
for p2i   = -k:k
for p3i  = -k:k
for p4i  = -k:k
                for p5i = -k:k
                    [p1, p2, p3, p4, p5, ds, du, da, dr] = deal(2^p1i, 2^p2i, 2^p3i, 2^p4i, 2^p5i,1,1,1,1);
                    sol = solve(subs(odes),'Real',true);
                    ss = [ss ; [p1i,p2i,p3i,p4i,p5i,length(sol.s) - 1]];
                    if length(sol.s) > 1
                        avg = [avg ; [p1i,p2i,p3i,p4i,p5i,mean(sol.u(2:length(sol.u)) + sol.a(2:length(sol.a)))]];
                    end;      end
            end
        end
    end
end



writematrix(double(avg),'avg.csv')
writematrix(double(ss), 'ss.csv')
