function [DMarray] = DMactuators2array(DMactuators)

DMarray = zeros(9,9);
DMarray(1,3:7) = DMactuators(1:5);
DMarray(2,2:8) = DMactuators(6:12);
% DMarray(3:7,1:9) = DMactuators(13:57);
for i=3:7
    DMarray(i,1:9) = DMactuators(13+(i-3)*9:13+(i-3)*9+8);
end
DMarray(8,2:8) = DMactuators(58:64);
DMarray(9,3:7) = DMactuators(65:69);