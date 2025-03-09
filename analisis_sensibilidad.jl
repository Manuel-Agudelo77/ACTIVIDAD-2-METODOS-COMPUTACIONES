using LinearAlgebra
using DataFrames
using CSV
using Plots
using SparseArrays

# Se realiza la ubicación de los CSV lines y nodes.

lines = DataFrame(CSV.File("actividad_2_Factor_de_sensibilidad/lines.csv"))
nodes = DataFrame(CSV.File("actividad_2_Factor_de_sensibilidad/nodes.csv"))

# Primero se calcula la matriz Bbus

function Matriz_Bbus(lines, nodes)
    """
    Calcula la matriz Bbus del flujo de carga DC.
    
    Entradas:   
        - lines: DataFrame con columnas FROM, TO, X.
        - nodes: DataFrame con columnas NUMBER, TYPE.
    
    Salidas:   
        - Bbus: Matriz de admitancias nodales sin el nodo slack.
    """
    total_nodes = nrow(nodes)
    total_lines = nrow(lines)
    Bbus = zeros(Float64, total_nodes, total_nodes)
    for k = 1:total_lines
        # Nodo de envío
        n1 = lines.FROM[k]
        # Nodo de recibo
        n2 = lines.TO[k]
        # Admitancia de la línea
        BL = 1/(lines.X[k])
        # Evitar errores por división entre cero
        if isinf(BL) || isnan(BL)
            error("Reactancia X en la línea $k es cero o inválida.")
        end
        @inbounds begin # Evita comprobaciones innecesarias en índices de Bbus.
        Bbus[n1,n1] += BL    # Dentro de la diagonal
        Bbus[n1,n2] -= BL    # Fuera de la diagonal
        Bbus[n2,n1] -= BL    # Fuera de la diagonal
        Bbus[n2,n2] += BL    # Dentro de la diagonal
        end
    end
    # Debo retirar el nodo slack
    slack = nodes[nodes.TYPE .== 3, "NUMBER"]
    if !isempty(slack)
        Bbus = Bbus[setdiff(1:end, slack), setdiff(1:end, slack)]
    end
    return Bbus
end


## Factores de sensibilidad lineales

# Calculo del Generation Shift Factor
function Generation_shift_α(lines, nodes, Bbus)
    """
    Calcula el Generation Shift Factor (GSF).
    
    Entradas:   
        - lines: DataFrame con columnas FROM, TO, X.
        - nodes: DataFrame con columnas NUMBER, TYPE.
        - Bbus: Matriz de susceptancias nodales (sin slack).

    Salidas:   
        - alpha_α: Matriz con los factores de sensibilidad por cambio de generación.
    """
    total_lines = nrow(lines)
    total_nodes = nrow(nodes)
    # Verificar nodo slack
    slack = nodes[nodes.TYPE .== 3, "NUMBER"]
    if length(slack) != 1
        error("Debe haber exactamente un nodo slack.")
    end
    s = slack[1]
    # Intentar invertir Bbus o usar la pseudo-inversa si es singular
    Bbus_inv = try
        Bbus_inv = inv(Bbus)
    catch
        Bbus_inv = pinv(Bbus)
    end
    # Expandir la matriz inversa para incluir el nodo slack
    row, col = size(Bbus_inv)
    W = zeros(row+1, col+1)

    # Añadiendo las filas y columnas del nodo slack    
    W[1:s-1, 1:s-1] = Bbus_inv[1:s-1, 1:s-1]  # Parte superior izquierda
    W[1:s-1, s+1:end] = Bbus_inv[1:s-1, s:end]  # Parte superior derecha
    W[s+1:end, 1:s-1] = Bbus_inv[s:end, 1:s-1]  # Parte inferior izquierda
    W[s+1:end, s+1:end] = Bbus_inv[s:end, s:end]  # Parte inferior derecha

    # Calculo del factor Generation Shift Factor
    alpha_α = zeros(total_lines, total_nodes)
    for i = 1:total_lines
        k = lines.FROM[i]
        m = lines.TO[i]
        x_i = lines.X[i]
        for j = 1:total_nodes
            alpha_α[i,j] = 1/(x_i)*((W[k,j]-W[m,j]))
        end
    end
    return alpha_α
end



# Cálculo del Factor de Distribución de Corte de Línea por Nodo (NBLODF)
function Node_Based_Line_Outage_Distribution_Factor(lines, nodes, Bbus)
    """
    Calcula el Factor de Distribución de Corte de Línea por Nodo (NBLODF).
    
    Entradas:   
        - lines: DataFrame con columnas FROM, TO, X.
        - nodes: DataFrame con columnas NUMBER, TYPE.
        - Bbus: Matriz de susceptancias nodales (sin slack).

    Salida:   
        - NBLODF: Matriz con los factores de sensibilidad por nodo ante la pérdida de una línea.
    """
    total_lines = nrow(lines)
    total_nodes = nrow(nodes)

    # Intentar invertir Bbus o usar la pseudo-inversa si es singular
    Bbus_inv = try
        inv(Bbus)
    catch
        pinv(Bbus)
    end

    # Obtener el nodo slack
    slack = nodes[nodes.TYPE .== 3, "NUMBER"]
    if length(slack) != 1
        error("Debe haber exactamente un nodo slack.")
    end
    s = slack[1]

    # Expandir la matriz inversa para incluir el nodo slack
    row, col = size(Bbus_inv)
    W = zeros(row+1, col+1)

    @views begin
        W[1:s-1, 1:s-1]   .= Bbus_inv[1:s-1, 1:s-1]  # Parte superior izquierda
        W[1:s-1, s+1:end] .= Bbus_inv[1:s-1, s:end]  # Parte superior derecha
        W[s+1:end, 1:s-1] .= Bbus_inv[s:end, 1:s-1]  # Parte inferior izquierda
        W[s+1:end, s+1:end] .= Bbus_inv[s:end, s:end] # Parte inferior derecha
    end

    # Cálculo del Factor de Distribución de Corte de Línea por Nodo (NBLODF)
    NBLODF = zeros(total_nodes, total_lines)
    for i = 1:total_lines
        k, m, x_i = lines.FROM[i], lines.TO[i], lines.X[i]
        den = x_i - (W[k, k] + W[m, m] - 2 * W[m, k])

        for j = 1:total_nodes
            NBLODF[j, i] = (x_i * (W[j, k] - W[j, m])) / den
        end
    end
    return NBLODF
end


# Calculo de los factores de distribución del corte de la linea
function Line_Outage_Distribution_Factor(lines, nodes, Bbus)
    """
    Calcula el Factor de Distribución de Corte de Línea (LODF).

    Entradas:
        - lines: DataFrame con columnas FROM, TO, X (impedancia).
        - nodes: DataFrame con columnas NUMBER, TYPE.
        - Bbus: Matriz de susceptancias nodales (sin slack).

    Salida:
        - beta: Matriz (total_lines × total_lines) con factores de distribución de corte de línea.
            Cada elemento beta[l, h] indica la sensibilidad del flujo en la línea `l`
            ante la desconexión de la línea `h`.

    Nota:
        - Se asume que hay un solo nodo slack en el sistema.
    """
    total_lines = nrow(lines)
    total_nodes = nrow(nodes)

    # Intentar invertir Bbus o usar la pseudo-inversa si es singular
    Bbus_inv = try
        inv(Bbus)
    catch
        pinv(Bbus)
    end
    # Expandir la matriz inversa para incluir el nodo slack
    row, col = size(Bbus_inv)
    W = zeros(row+1, col+1)
    # Obtener el nodo slack
    slack = nodes[nodes.TYPE .== 3, "NUMBER"]
    if length(slack) != 1
        error("Error: Se encontraron $(length(slack)) nodos slack, pero debe haber exactamente uno.")
    end
    s = slack[1]
    @views begin
    W[1:s-1, 1:s-1]   .= Bbus_inv[1:s-1, 1:s-1]  # Parte superior izquierda
    W[1:s-1, s+1:end] .= Bbus_inv[1:s-1, s:end]  # Parte superior derecha
    W[s+1:end, 1:s-1] .= Bbus_inv[s:end, 1:s-1]  # Parte inferior izquierda
    W[s+1:end, s+1:end] .= Bbus_inv[s:end, s:end] # Parte inferior derecha
    end
    # Factores de distribución del corte de la linea
    beta = Matrix{Float64}(undef, total_lines, total_lines)
    for l = 1:total_lines
        k = lines.FROM[l]
        m = lines.TO[l]
        x_i = lines.X[l]
        for h = 1:total_lines
            i = lines.FROM[h]
            j = lines.TO[h]
            x_j = lines.X[h]
            if l != h
                beta[l,h] = (x_j/x_i)*(W[i,k]-W[j,k]-W[i,m]+W[j,m])/(x_j-W[i,i]-W[j,j]+2*W[i,j])
            end
        end
    end
    return beta
end

Bbus = Matriz_Bbus(lines, nodes)
display(Line_Outage_Distribution_Factor(lines, nodes, Bbus))