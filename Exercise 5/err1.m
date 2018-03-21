%error analysis script for relaxation2

close all;
rfinal = 0.01;
rstep = 0.00025;

for rmax = 0:rstep:rfinal
    relaxation2a;
    error = max(max(err));
    plot(log(rmax),log(error),'r+');
    hold on;

end

xlabel("log(maximum value of residue)");
ylabel("log(maximum error)");