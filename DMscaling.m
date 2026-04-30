function [DMvector_scl, DMarray_scl_fill] = DMscaling(DMvector)

DMarray = DMactuators2array(DMvector);

DMarray_scl = imresize(DMarray,0.5);
DMarray0 = zeros(9,9);
DMarray_scl_fill = DMarray0;
DMarray_scl_fill(3:7,3:7) = DMarray_scl;
DMarray_scl_fill(DMarray_scl_fill==0 )= nan;
DMarray_scl_fill=fillmissing(DMarray_scl_fill,'nearest',1);
DMarray_scl_fill=fillmissing(DMarray_scl_fill,'nearest',2);

DMvector_scl = DMarray2actuators(DMarray_scl_fill)';

