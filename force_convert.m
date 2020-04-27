function [ff]=force_convert(fsrc,thetas,phis)
ff(1)=sin(thetas)*cos(phis)*fsrc(1)+sin(thetas)*sin(phis)*fsrc(2)+cos(thetas)*fsrc(3);
ff(2)=cos(thetas)*cos(phis)*fsrc(1)+cos(thetas)*sin(phis)*fsrc(2)-sin(thetas)*fsrc(3);
ff(3)=-sin(phis)*fsrc(1)+cos(phis)*fsrc(2);
%ff=cart2sphvec(fsrc,thetas,phis)
end