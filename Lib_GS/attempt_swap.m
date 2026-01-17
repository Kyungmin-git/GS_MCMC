function [A,B] = attempt_swap(A,B)
    if ~isfinite(A.ll_amp) || ~isfinite(B.ll_amp)
        return;
    end

    llA = A.ll_amp + A.ll_ph;
    llB = B.ll_amp + B.ll_ph;

    logalpha = (1/A.Temp - 1/B.Temp) * (llB - llA);
    logalpha = min(0, logalpha);

    if log(rand()) < logalpha
        payloadA = get_payload(A);
        payloadB = get_payload(B);

        A = set_payload(A, payloadB);
        B = set_payload(B, payloadA);
    end
end

function P = get_payload(C)
    P.m = C.m;
    P.state = C.state;
    P.model = C.model;
    P.ll_amp = C.ll_amp;
    P.ll_ph  = C.ll_ph;
    P.lp     = C.lp;
end

function C = set_payload(C, P)
    C.m = P.m;
    C.state = P.state;
    C.model = P.model;
    C.ll_amp = P.ll_amp;
    C.ll_ph  = P.ll_ph;
    C.lp     = P.lp;
end
