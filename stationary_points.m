syms s u a r p1 p2 p3 p4 p5 ds du da dr


ode1 =  (-s*p1*(u+a)+u*p2*(s+r)  -ds*s + 1)-s *(1 - ds*s - du*u - da*a - dr*r);
ode2 =  (s*p1*(u+a)-u*p2*(s+r) - u*p3 -du*u )-u*(1 - ds*s - du*u - da*a - dr*r) ;
ode3 =  (u*p3 -a*p4*(s+r) + r*p5*(u+a) -da*a )-a*(1 - ds*s - du*u - da*a - dr*r);
ode4 =  (a*p4*(s+r) - r*p5*(u+a) -dr*r )-r*(1 - ds*s -du*u - da*a - dr*r);
sums = s+u+a+r;


odes = [ode1 == 0, ode2 == 0, ode3 == 0 , ode4 == 0, sums  == 1];

o = [ode1, ode2, ode3 , ode4];

j = jacobian(o,[s,u,a,r]);


h = 0.1;

res = [];


for s = 0:h:1
    for u = 0:h:1
        for a = 0:h:1
            r = 1 - s - u - a;
            if r >= 0
                    
                    x = [s;u;a;r];
                    d = (subs(j) * x +x);
                    v = [s,u,a,r,d(1),d(2),d(3),d(4)];
                    res = [res;v];
            end
        end
    end
end

clear s u a r

[p1, p2, p3, p4, p5, ds, du, da, dr] = deal(10,1,1, 1, 1,0.1,0.5,1,0.1);

sres = subs(res) ;

sol = solve(subs(odes),'Real',true);


[s,u,a,r] = deal(0.122386,0.495595,0.339047,0.042971);
chnage = simplify(subs(o));

[s,u,a,r] = deal(1,0,0,0);
real(eig(subs(j)))


