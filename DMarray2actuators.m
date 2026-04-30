function [DMactuators] = DMarray2actuators(DMarray)

DMactuators = zeros(69,1);

DMactuators(1:5) = DMarray(1,3:7);
DMactuators(6:12) = DMarray(2,2:8);
% DMactuators(13:57) = DMarray(3:7,1:9);
for i=3:7
    DMactuators(13+(i-3)*9:13+(i-3)*9+8) = DMarray(i,1:9);
end
DMactuators(58:64) = DMarray(8,2:8);
DMactuators(65:69) = DMarray(9,3:7);
