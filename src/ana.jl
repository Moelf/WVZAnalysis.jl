const Z_m = 91.1876 * 10^3 # in MeV

mass(lv::LorentzVector) = sqrt(dot(lv, lv))

function Find_Z_Pairs(v_l_pids, v_l_tlv, v_l_wgt)
    v_Z_pair = Tuple{Int,Int}[]
    v_Z_wgt = Float64[]

    v_ignore = Set{Int}()
    @inbounds while length(v_ignore) < length(v_l_pids)
        M = 0.0
        local temp_tup
        for i in eachindex(v_l_pids) # electrons loop
            i ∈ v_ignore && continue
            for j in (i + 1):length(v_l_pids)
                j ∈ v_ignore && continue
                v_l_pids[i] != -v_l_pids[j] && continue # require OS
                m0 = mass(v_l_tlv[i] + v_l_tlv[j])
                if abs(m0 - Z_m) < abs(M - Z_m)
                    M = m0
                    temp_tup = (i, j)
                end
            end
        end
        iszero(M) && break # can't find any more pairs
        push!(v_ignore, temp_tup[1])
        push!(v_ignore, temp_tup[2])

        push!(v_Z_pair, temp_tup)
        push!(v_Z_wgt, v_l_wgt[temp_tup[1]] * v_l_wgt[temp_tup[2]])
    end

    return v_Z_pair, v_Z_wgt, v_ignore
end

function Find_m4l(v_Z_pair, v_l_tlv, v_l_order)
    m4l = first(v_Z_pair)

    @inbounds for i in eachindex(v_l_order)
        length(m4l) >= 4 && break
        (v_l_order[i] == m4l[1] || v_l_order[i] == m4l[2]) && continue
        m4l = (m4l..., v_l_order[i])
    end
    # require 2SFOS
    length(m4l) != 4 && return -1.0

    tlv_4l = zero(LorentzVector)
    for idx in m4l
        tlv_4l += v_l_tlv[idx]
    end

    return mass(tlv_4l)
end

@pour initial_cut begin
    v_l_pid = vcat(evt.v_e_pid, evt.v_m_pid)
    nlepton = length(v_l_pid)
    nlepton <= 3 && return false

    v_l_tlv = vcat(evt.v_e_tlv, evt.v_m_tlv)
    v_l_wgt = vcat(evt.v_e_wgtLoose, evt.v_m_wgtLoose)

    v_Z_pair, v_Z_wgt, v_ignore = Find_Z_Pairs(v_l_pid, v_l_tlv, v_l_wgt)
    length(v_Z_wgt) == 0 && return false

    best_p = first(v_Z_pair)
    abs(mass(v_l_tlv[best_p[1]] + v_l_tlv[best_p[2]]) - Z_m) > 20e3 && return false

    v_l_order = sortperm(v_l_tlv; by=mass)
    mass_4l = Find_m4l(v_Z_pair, v_l_tlv, v_l_order)
end # end of initial_cut

function main_looper(mytree)
    for evt in mytree
        ### initial_cut
        pass_initial_cut = @initial_cut
        ### ZZZ_cut
        pass_ZZZ_cut = @ZZZ_cut

        break
    end
end
