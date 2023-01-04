classdef BoundaryOptimizedSBPOps
    % Implements the boundary-optimized SBP operators for first derivatives (D1)
    % and second derivative variable coefficient (D2)
    %
    % The boundary closure of D2 uses the first and last rows of D1
    % i.e. the operators are fully compatible.
    %
    % The D1 operators are presented in
    %   Mattsson, Ken and Almquist, Martin and van der Weide, Edwin
    %   Boundary optimized diagonal-norm {SBP} operators
    %   Journal of Computational Physics 2018
    %   https://doi.org/10.1016/j.jcp.2018.06.010
    % 
    % 
    properties
        H % Norm (quadrature) matrix
        HI % H^-1
        D1 % SBP operator approximating first derivative d/dx
        D2 % SBP operator approximating second derivative d/dx(c)d/dx
        DI % Dissipation operator
        e_l % Left boundary operator
        e_r % Right boundary operator
        d1_l % Left boundary first derivative
        d1_r % Right boundary first derivative
        m % Number of grid points.
        h % Step size
        x % grid
        options % Struct holding options used to create the operators
    end

    methods

        function obj = BoundaryOptimizedSBPOps(m, lim, order, options)
            % BoundaryOptimizedSBPOps BoundaryOptimizedSBPOps(m, lim, order, options)
            %   m - number of gridpoints
            %   lim - cell array holding the limits of the domain
            %   order - order of the operator
            %   options - struct holding options used to construct the operators
            %             options.stencil_width:   {'narrow', 'seminarrow', 'wide'}
            %                 Specifies the interior stencil width of D2.
            %                 narrow: narrow-stencil (default)
            %                 seminarrow: two additional stencil points compared to narrow
            %                 wide: wide-stencil obtained by applying D1 twice
            %             options.variable_coeffs: {true, false}
            %                 true: obj.D2 is a function handle D2(c) returning a matrix 
            %                       for coefficient vector c (default)
            %                 false: obj.D2 is a matrix.
            %             options.AD: {'op', 'upwind'}
            %                  Specifies how the artificial dissipation operator DI is constructed
            %                 'op': order-preserving AD (preserving interior stencil order) (default)
            %                 'upwind': upwind AD (order-1 upwind interior stencil)
            util.default_arg('options', struct);
            util.default_field(options,'stencil_width','narrow');
            util.default_field(options,'variable_coeffs',true);
            util.default_field(options,'AD','op');
            assert(any(order == 4:2:12),'order must be one of 4, 6, 8, 10, 12');
            [x, h] = implementations.boundaryOptimizedGrid(lim, m, order);
            switch order
                case 4
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = implementations.bopt_ops_4(m, h, options);
                case 6
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = implementations.bopt_ops_6(m, h, options);
                case 8
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = implementations.bopt_ops_8(m, h, options);
                case 10
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = implementations.bopt_ops_10(m, h, options);
                case 12
                    [obj.H, obj.HI, obj.D1, obj.D2, obj.DI] = implementations.bopt_ops_12(m, h, options);
            end

            if ~options.variable_coeffs
                obj.D2 = obj.D2(ones(m,1));
            end

            % Boundary operators
            obj.e_l = sparse(m, 1); obj.e_l(1) = 1;
            obj.e_r = sparse(m, 1); obj.e_r(m) = 1;
            obj.d1_l = (obj.e_l' * obj.D1)';
            obj.d1_r = (obj.e_r' * obj.D1)';
            % grid data
            obj.x = x;
            obj.h = h;
            obj.m = m;

            % misc
            obj.options = options;
        end
    end

end
