function grd = orbitTransferGrd(Z)
% computes the gradient

output = orbitTransferObj_Jac(Z);
grd    = output;

end

