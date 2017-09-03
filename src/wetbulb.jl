


function calc_W_from_B(Tk, B, P, EPS=1e-8, MAXITER=100)

    pws = Pws(B)
    efac = efactor(B,P)
    xsv = efac * pws / P
    w2 = Mv/Ma * xsv/(1-xsv)


    w = (enthalpyair(B,P) - enthalpyair(Tk,P) - w2*(enthalpywi(B) - enthalpyvapor(B))) /
        (enthalpyvapor(Tk) - enthalpywi(B))

    for iter = 1:MAXITER
        f = aux_WB(w, Tk, B, P, pws, efac)
        df = (aux_WB(w+1e-4*w2, Tk, B, P, pws, efac) - f) / (1e-4*w2)
        dw = -f/df
        w = w + dw
        println(dw)

        if abs(dw) < EPS*w2
            return w
        end

    end

    return w
end

function aux_WB(w, Tk, B, P, pws, efac)
    xv1 = w / (Mv/Ma+w)
    xv2 = efac * pws / P
    w2 = Mv / Ma * xv2 / (1.0-xv2)
    #(1.0+w)*enthalpymoist(Tk,P,xv1) + (w2-w)*enthalpywi(B) - (1.0+w2)*enthalpymoist(B,P,xv2)
    enthalpymoist(Tk,P,xv1) + (w2-w)*enthalpywi(B) - enthalpymoist(B,P,xv2)
end

function aux_WB(w, Tk, B, P)
    pws = Pws(B)

    efac = efactor(B,P)
    
    xv1 = w / (Mv/Ma+w)
    xv2 = efac * pws / P
    w2 = Mv / Ma * xv2 / (1.0-xv2)
    #(1.0+w)*enthalpymoist(Tk,P,xv1) + (w2-w)*enthalpywi(B) - (1.0+w2)*enthalpymoist(B,P,xv2)
    enthalpymoist(Tk,P,xv1) + (w2-w)*enthalpywi(B) - enthalpymoist(B,P,xv2)
end

function calcwetbulb(Tk, P, xv, EPS=1e-8, MAXITER=200)

    w = humrat(xv)

    B = Tk - 1.0 # Initial guess
    h = 1e-7
    i = 0
    dB = 0.0
    for i = 1:MAXITER
        f = aux_WB(w, Tk, B, P)
        df = (aux_WB(w, Tk, B + h, P) - f) / h
        dB = -f/df
        B = B + dB

        if abs(dB) < EPS
            return B
        end
    end
    throw(ConvergenceError("Wet Bulb temperature failed to converge!", B, i, dB))
    return B  # Later on I will have to check covergence.
      
end
