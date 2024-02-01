function [Fsee, Fce, Fpee] = findMuscleForce(var, musvar, l_ce, v_ce, a, part_only)

    if nargin == 5 % Not checking force-length or force-velocity relationship
        for i = 1:length(l_ce)
            Fpee(i) = findFpee(l_ce(i), var, musvar);
            Fce(i) = findFce(l_ce(i), a(i), v_ce(i), var);
            
            Fsee(i)= Fpee(i) + Fce(i);
        end
    else %check force-length or force-velocity relationship
        if strcmpi(part_only, 'Fpee')
            for i = 1:length(l_ce)
                Fsee(i) = a(i)*findFpee(l_ce(i), var, musvar);
            end
        elseif strcmpi(part_only, 'flce')
            for i = 1:length(l_ce)
                Fsee(i) = a(i)*find_flce(l_ce(i), var);
            end
        elseif strcmpi(part_only, 'flcePEE')
            for i = 1:length(l_ce)
                Flce = find_flce(l_ce(i), var);
                Fpee = findFpee(l_ce(i), var, musvar);
                Fsee(i) = a(i)*(Flce+Fpee);
            end
        elseif strcmpi(part_only, 'gvce')
            for i = 1:length(v_ce)
                Fsee(i) = a(i)*find_gvce(v_ce(i), var);
            end
        end
    end
end

function Fpee = findFpee(lce, var, musvar)
    Fpee = 0;
    
    xPE = lce - var.PEEslack;
    if xPE > 0
        Fpee = Fpee + var.kPEE*(xPE*musvar.l_opt)^2;
    end
end

function Fce = findFce(lce, a, vce, var)
    %Returns the contractile element force, normalized to fmax
    F1 = find_flce(lce, var);
    F2 = find_gvce(vce, var);
    
    Fce = a*F1*F2;
end

function F1 = find_flce(lce, var)
    x = (lce - 1.0)/var.W; %width of the force-length relationship
    F1 = exp(-x*x);
end

function F2 = find_gvce(vce, var)
    if vce < 0
        F2 = (var.v_max + vce)/(var.v_max - vce/var.Arel);
    else
        F2 = (var.gmax*vce + var.c3)/(vce+var.c3);
    end
end