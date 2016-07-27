module MyRoots
#=
	Julia module to find roots of equations using some standard methods.
=#

function bisect(f,a,b; TOL=1e-9, MAXN=3000, del=1e-9)
    # computes f(x)=0 using bisection algorithm
    u = f(a);
    v = f(b);
    e = b-a;
    if sign(u)==sign(v)
        throw(DomainError());
    end
    for k=1:MAXN
        e = e/2;
        c = a+e;
        w = f(c);
        if (abs(e)<del) || (abs(w)<TOL)
            return c;
        end
        if sign(w)!=sign(u)
            b = c;
            v = w;
        else
            a = c;
            u = w;
        end
    end
end

function newton(p0,f,df; TOL=1e-9, MAXN=3000)
    # Function to calculate root f(x)=0 using Newton's method.
    # Requires initial guess p0 and functions f, df.
    # Optional tolerance and max-iter settings.
    # Returns root if found, throws bound error otherwise.
  for i=1:MAXN
    while i<=MAXN
      p=p0 - f(p0)/df(p0);
      if abs(p-p0)<TOL
        return p
      end
      i+=1;
      p0=p;
    end
  end
  throw(BoundsError);
end

function secant(p0,p1,f;TOL=1e-5, MAXN=2000)
    #= Func
    =#
    for i=1:MAXN
      while i<=MAXN
        p=p1 - f(p1)*(p1-p0)/(f(p1)-f(p0));
        if abs(p-p0)<TOL
          return p
        end
        i+=1;
        p0=p1;
        p1=p;
      end
    end
    throw(BoundsError);
end

end
