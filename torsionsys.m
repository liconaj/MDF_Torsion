classdef torsionsys
    %TORSIONSYS solucionador problema de torsión en una viga
    %   Solucionador de problema de torsión en vigas de sección regular

    properties (GetAccess = public)
        % Propiedades públicas, accesibles desde fuera de la clase
        mesh        % Malla
        meshindex   % Índice de la malla
        sketch      % Boceto
        matA        % Matriz A
        vecF        % Vector F
        Phi         % Función Phi
        imgPhi      % Imagen de Phi
        theta       % Ángulo de deformación por unidad de distancia
        gamma       % Ángulo de deformación total
        tauxz       % Componente tau_xz
        tauyz       % Componente tau_yz
        tau         % Magnitud de tau
        tauMax      % Máximo esfuerzo cortante
        ny          % Factor de seguridad
        J, W, L, G, h, T, sy  % Parámetros y propiedades del sistema
    end

    properties (GetAccess = private)
        % Propiedades privadas, solo accesibles dentro de la clase
        targetT         % Valor objetivo de T (Torque) 
        profile_data    % Datos del perfil
        n               % Número de nodos
        m               % Tamaño de matriz
    end

    methods
        function obj = torsionsys(profile_data, params)
            %TORSIONSYS Crear instancia de la clase
            %   Inicializa las propiedades con datos dados
            obj.profile_data = profile_data;
            obj.W = params.W;
            obj.L = params.L;
            obj.G = params.G;
            obj.h = params.h;
            obj.sy = params.sy;
            obj.targetT = params.T;
            obj = obj.setupemesh(); % Configura la malla
            obj.J = obj.calcJ();    % Calcula el momento polar de inercia
        end

        function obj = setupsystem(obj)
            %SETUPSYSTEM Configuración del sistema de ecuaciones: define la
            %matriz A para resolver el problema

            % Inicializa la matriz A con valores específicos
            obj.matA = -4*speye(obj.n);
            % Recorre la malla y asigna coeficientes a la matriz A
            for y = 1:obj.m
                for x = 1:obj.m
                    ii = obj.meshindex(y, x);
                    if isnan(ii) || ii == 0
                        continue
                    end
                    ncase = obj.mesh(y-1:y+1, x-1:x+1);
                    adjcoefs = obj.getadjcoefs(ncase);
                    jj(1) = obj.meshindex(y, x-1);
                    jj(2) = obj.meshindex(y-1, x);
                    jj(3) = obj.meshindex(y, x+1);
                    jj(4) = obj.meshindex(y+1, x);
                    for idx = 1:length(adjcoefs)
                        c = adjcoefs(idx);
                        if c > 0
                            obj.matA(ii, jj(idx)) = c;
                        end
                    end
                end
            end
        end

        function obj = solvesystem(obj)
            %SOLVESYSTEM Resuelve el sistema de ecuaciones para encontrar
            %la distribución de Phi y otros valores.
            eq = @(T) obj.calcT(T) - obj.targetT;
            Ttheta = fsolve(eq, obj.targetT, optimoptions(@fsolve, "Display", "off"));
            obj.theta = Ttheta / (obj.J * obj.G);
            obj.gamma = obj.theta * obj.L;
            obj.Phi = obj.calcPhi(obj.theta);
            obj.T = 2*sum(obj.Phi) * obj.h^2;
        end

        function obj = arrangesolution(obj)
            %ARRANGESOLUTION Organiza la distribución de Phi en la malla
            obj.imgPhi = obj.mesh;
            for y = 1:obj.m
                for x = 1:obj.m
                    idx = obj.meshindex(y, x);
                    if isnan(idx) || idx == 0
                        continue
                    end
                    obj.imgPhi(y, x) = obj.Phi(idx);
                end
            end
            % invertir columnas para visualización
            obj.imgPhi = obj.imgPhi(:,end:-1:1);
            % Calcular las componentes tau_xz, tau_yz
            obj.tauxz = nan(obj.m-2);
            obj.tauyz = nan(obj.m-2);
            for y = 2:obj.m-1
                for x = 2:obj.m-1
                    if isnan(obj.imgPhi(y,x))
                        continue
                    end
                    dPhix = obj.imgPhi(y,x+1)-obj.imgPhi(y,x-1);
                    obj.tauxz(y-1,x-1) = dPhix/(2*obj.h);
                    dPhiy = obj.imgPhi(y+1,x)-obj.imgPhi(y-1,x);
                    obj.tauyz(y-1,x-1) = -dPhiy/(2*obj.h);
                end
            end
            % Calcula las magnitudes de tau
            obj.tau = (obj.tauxz.^2+obj.tauyz.^2).^0.5;
            obj.tauMax = max(obj.tau, [], "all");
            obj.ny = obj.sy / (2*obj.tauMax);
            if obj.tauMax > obj.sy/2
                warning("Sistema falla estáticamente por criterio de Tresca")
            end
        end

        function Phi = calcPhi(obj, theta)
            %CALCPHI Calcular el vector solución Phi
            F = -2*obj.G*theta;
            obj.vecF = F * ones(obj.n, 1);
            Phi = obj.matA\obj.vecF;
        end

        function T = calcT(obj, Ttheta)
            %CALCT Calcular el torque T a partir del vector Phi dada por
            %Ttheta
            theta_guess = Ttheta / (obj.J * obj.G);
            T = 2*sum(obj.calcPhi(theta_guess)) * obj.h^2;
        end

        function showPhi(obj)
            %SHOWPHI Mostrar distribución de Phi en la geometría
            figure
            xt = (1:obj.m)*obj.h;
            yt = (1:obj.m)*obj.h;
            contourf(xt,yt,obj.imgPhi, 20, 'LineColor', 'flat')
            title("Distribución \Phi")
            colorbar
            daspect([1 1 1]) % Ajustar proporción ejes de la gráfica
            % set(gca,'Color','k') % Cambia el color de fondo (opcional)
        end

        function showtau(obj)
            %SHOWTAU Mostrar magnitud esfuerzos cortantes en geometría
            figure
            xt = (2:obj.m-1)*obj.h;
            yt = (2:obj.m-1)*obj.h;
            contourf(yt,xt,obj.tau, 20, 'LineColor', 'flat')
            title("Distribución magnitud \tau")
            hc = colorbar;
            title(hc, "MPa") % Mostrar unidad esfuerzo
            colormap turbo
            daspect([1 1 1]) % Ajustar proporción ejes de la gráfica
        end

        function showtaucamp(obj)
            %SHOWTAUCAMP Mostrar vectores de esfuerzos cortantes en la
            %geometría
            figure
            p = (obj.m-3)/17;
            xt = (2:p:obj.m-1)*obj.h;
            yt = (2:p:obj.m-1)*obj.h;
            U = obj.tauyz(floor(1:p:end),floor(1:p:end));
            V = obj.tauxz(floor(1:p:end),floor(1:p:end));
            quiver(xt,yt,U,V)
            hold on
            xt = (2:obj.m-1)*obj.h;
            yt = (2:obj.m-1)*obj.h;
            contour(xt,yt,obj.tau,20, 'LineColor', 'flat')
            title("Distribución vectores \tau")
            hc = colorbar;
            title(hc, "MPa") % Mostrar unidad esfuerzo
            colormap turbo
            daspect([1 1 1]) % Ajustar proporción ejes de la gráfica
        end
        function showtaucomps(obj)
            %SHOWTAUCOMPS Mostrar magnitudes de las componentes de los
            %esfuerzos cortantes
            figure
            xt = (2:obj.m-1)*obj.h;
            yt = (2:obj.m-1)*obj.h;
            
            subplot(1,2,1)
            contourf(xt,yt,obj.tauxz(end:-1:1,:), 20, 'LineColor', 'flat')
            title("\tau_xz")
            hc = colorbar;
            title(hc, "MPa") % Mostrar unidad esfuerzo
            colormap turbo
            daspect([1 1 1]) % Ajustar proporción ejes de la gráfica
            
            subplot(1,2,2)
            contourf(xt,yt,obj.tauyz(end:-1:1,:), 20, 'LineColor', 'flat')
            title("\tau_yz")
            hc = colorbar;
            title(hc, "MPa") % Mostrar unidad esfuerzo
            colormap turbo
            daspect([1 1 1]) % Ajustar proporción ejes de la gráfica
        end
    end
    methods (Access = private)
        function obj = setupemesh(obj)
            % Función privada para configurar la malla de puntos
            size = length(obj.profile_data);
            cellsize = obj.W/size; % Tamaño de división
            scale = cellsize/obj.h; % Escala para el tamaño de la malla
            % Escalar datos a número de celdas
            obj.sketch = imresize(obj.profile_data, scale, "nearest");
            obj.m = length(obj.sketch) + 1;
            rawmesh = nan(obj.m); % Inicialización de la malla
            % Creación de la malla a partir del celdas
            for y = 1:obj.m
                for x = 1:obj.m
                    if x == 1
                        xs = x;
                    elseif x == obj.m
                        xs = x - 1;
                    else
                        xs = [x-1, x];
                    end
                    if y == 1
                        ys = y;
                    elseif y == obj.m
                        ys = y - 1;
                    else
                        ys = [y-1, y];
                    end
                    neighbors = obj.sketch(ys, xs); % Vecinos en el perfil
                    val = sum(neighbors, "all"); % Suma de los vecinos
                    if val > 2 % && ~any(isnan(neighbors(:)))
                        rawmesh(y,x) = 1;
                    elseif val > 0
                        rawmesh(y,x) = 0;
                    end
                end
            end
            obj.mesh = rawmesh;
            obj.meshindex = zeros(obj.m);
            obj.n = 0;
            idx = 0;
            for y = 1:obj.m
                for x = 1:obj.m
                    if ~isequal(rawmesh(y,x), 1)
                        continue
                    end
                    idx = idx + 1;
                    obj.n = obj.n + 1;
                    sumy = sum(rawmesh(y,[x-1, x+1]));
                    sumx = sum(rawmesh([y-1, y+1], x));
                    obj.mesh(y,x) = sumy + sumx;
                    obj.meshindex(y,x) = idx;
                end
            end
        end

        function J = calcJ(obj)
            % Función privada para calcular el momento polar de inercia J
            J = 0; % Inicialización del momento polar de inercia
            dA = obj.h^2; % Área diferencial dA
            cx = obj.W/2; % Coordenada x del centro del perfil
            cy = obj.W/2; % Coordenada y del centro del perfil
            [sm, sn] = size(obj.sketch); % Número de dA discretizadas
            for jj = 1:sm
                for ii = 1:sn
                    if ~obj.sketch(jj,ii)
                        continue
                    end
                    x = ii*obj.h; % Coordenada x de dA en la malla
                    y = jj*obj.h; % Coordenada y de dA en la malla
                    % Distancia al cuadrado entre dA y el centro del perfil
                    r2 = (x-cx)^2+(y-cy)^2;
                    J = J + r2 * dA;
                end
            end
        end
    end
    methods (Static)
        function coefs = getadjcoefs(ncase)
            % Obtener coeficientes de nodos adyacentes
            ntype = ncase(2,2);  % Tipo de nodo
            switch ntype
                case 4
                    coefs = [1 1 1 1];
                case 3
                    coefs = torsionsys.shiftcoefs(ncase, [0 1 1 1], 0, 0);
                case 2
                    coefs = torsionsys.shiftcoefs(ncase, [1 1 0 0], 4, 1);
            end
        end
        function coefs = shiftcoefs(nodescase, defcoefs, nodeseek, seekcorner)
            % Cambiar coeficientes según rotación del caso
            % Explicación:
            % adjnodes => [1 2 3 4]
            % "side"  . 2 .   "corner"   1 . 2
            %         1 . 3              . . .
            %         . 4 .              4 . 3
            adjnodes = zeros(1, 4);
            y = 2; x = 2; % Coordenadas nodo central en caso
            if seekcorner
                adjnodes(1) = nodescase(y-1, x-1); % Nodo superior izquierdo
                adjnodes(2) = nodescase(y-1, x+1); % Nodo superior derecho
                adjnodes(3) = nodescase(y+1, x+1); % Nodo inferior derecho
                adjnodes(4) = nodescase(y+1, x-1); % Nodo inferior izquierdo
            else
                adjnodes(1) = nodescase(y, x-1); % Nodo izquierdo
                adjnodes(2) = nodescase(y-1, x); % Nodo superior
                adjnodes(3) = nodescase(y, x+1); % Nodo derecho
                adjnodes(4) = nodescase(y+1, x); % Nodo inferior
            end
            shiftsteps = find(adjnodes == nodeseek) - 1;
            coefs = circshift(defcoefs, shiftsteps);
        end
    end
end