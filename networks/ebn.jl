# using Distributed
# numProcs = 10
# addprocs(numProcs)

# @everywhere begin
using EnhancedBayesianNetworks
using EnhancedBayesianNetworks: evaluate!
using CSV
using DataFrames
using JLD2
using Dates

const SIMULATIONS = 1 # number of Monte Carlo simulations
const threshold = 1243.9

const current_dir = pwd()

``` PGA node peak ground acceleration
```
cpt_pga = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:PGA)
pga_probs = CSV.read("networks/csv/PGA_probability.csv", DataFrame)[:, 1]
pga_states = ["PGA_$(lpad(i, 2, '0'))" for i in 0:length(pga_probs)]
map((p, st) -> cpt_pga[:PGA=>Symbol(st)] = p, pga_probs, pga_states)
pga_node = DiscreteNode(:PGA, cpt_pga)


``` AGE node
```
cpt_age = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(:AGE)
age_states = ["AGE_$(i)0" for i in 0:5]
map(st -> cpt_age[:AGE=>Symbol(st)] = 1 / length(age_states), age_states)
age_node = DiscreteNode(:AGE, cpt_age)


``` LOCA node
```
data = CSV.read("networks/csv/LOCA_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
rename!(data, :AGE_0 => :AGE_00)
df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
df_long2 = deepcopy(df_long1)
df_long2[!, :Π] = 1 .- df_long2[!, :Π]
df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
insertcols!(df_long1, 3, :LOCA => fill(:YES_LOCA, nrow(df_long1)))
insertcols!(df_long2, 3, :LOCA => fill(:NO_LOCA, nrow(df_long2)))
df_loca = sort(vcat(df_long1, df_long2))
loca_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_loca)
loca_node = DiscreteNode(:LOCA, loca_cpt)


``` LOCA-TIME node
```
t_loca_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOCA)
t_loca_cpt[:LOCA=>:YES_LOCA] = LogNormal(3.3, 1)
t_loca_cpt[:LOCA=>:NO_LOCA] = Normal(1200, 0)
t_loca_node = ContinuousNode(:t_loca, t_loca_cpt)


``` LOOP node
```
data = CSV.read("networks/csv/LOOP_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
df_loop = sort(vcat(
    DataFrame(PGA=Symbol.(data.PGA), LOOP=:YES_LOOP, Π=Real.(data.LOOP)),
    DataFrame(PGA=Symbol.(data.PGA), LOOP=:NO_LOOP, Π=Real.(1 .- data.LOOP))
))
loop_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_loop)
loop_node = DiscreteNode(:LOOP, loop_cpt)


``` LOOP-TIME node
```
t_loop_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOOP)
t_loop_cpt[:LOOP=>:YES_LOOP] = LogNormal(3.3, 1)
t_loop_cpt[:LOOP=>:NO_LOOP] = Normal(1200, 0)
t_loop_node = ContinuousNode(:t_loop, t_loop_cpt)


``` MSLB node
```
data = CSV.read("networks/csv/MSLB_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
rename!(data, :AGE_0 => :AGE_00)
df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
df_long2 = deepcopy(df_long1)
df_long2[!, :Π] = 1 .- df_long2[!, :Π]
df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
insertcols!(df_long1, 3, :MSLB => fill(:YES_MSLB, nrow(df_long1)))
insertcols!(df_long2, 3, :MSLB => fill(:NO_MSLB, nrow(df_long2)))
df_mslb = sort(vcat(df_long1, df_long2))
mslb_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_mslb)
mslb_node = DiscreteNode(:MSLB, mslb_cpt)


``` MSLB-TIME node
```
t_mslb_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:MSLB)
t_mslb_cpt[:MSLB=>:YES_MSLB] = LogNormal(3.3, 1)
t_mslb_cpt[:MSLB=>:NO_MSLB] = Normal(1200, 0)
t_mslb_node = ContinuousNode(:t_mslb, t_mslb_cpt)


``` MSLBH2 node
```
MSLBH2_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}([:PGA, :MSLBH2])
MSLBH2_cpt[:PGA=>:PGA_00, :MSLBH2=>:YES_MSLBH2] = 0
MSLBH2_cpt[:PGA=>:PGA_00, :MSLBH2=>:NO_MSLBH2] = 1
MSLBH2_cpt[:PGA=>:PGA_01, :MSLBH2=>:YES_MSLBH2] = 0.045802235
MSLBH2_cpt[:PGA=>:PGA_01, :MSLBH2=>:NO_MSLBH2] = 1 - 0.045802235
MSLBH2_cpt[:PGA=>:PGA_02, :MSLBH2=>:YES_MSLBH2] = 0.290714844
MSLBH2_cpt[:PGA=>:PGA_02, :MSLBH2=>:NO_MSLBH2] = 1 - 0.290714844
MSLBH2_cpt[:PGA=>:PGA_03, :MSLBH2=>:YES_MSLBH2] = 0.5450015
MSLBH2_cpt[:PGA=>:PGA_03, :MSLBH2=>:NO_MSLBH2] = 1 - 0.5450015
MSLBH2_cpt[:PGA=>:PGA_04, :MSLBH2=>:YES_MSLBH2] = 0.72052405
MSLBH2_cpt[:PGA=>:PGA_04, :MSLBH2=>:NO_MSLBH2] = 1 - 0.72052405
MSLBH2_cpt[:PGA=>:PGA_05, :MSLBH2=>:YES_MSLBH2] = 0.828946831
MSLBH2_cpt[:PGA=>:PGA_05, :MSLBH2=>:NO_MSLBH2] = 1 - 0.828946831
MSLBH2_cpt[:PGA=>:PGA_06, :MSLBH2=>:YES_MSLBH2] = 0.894119709
MSLBH2_cpt[:PGA=>:PGA_06, :MSLBH2=>:NO_MSLBH2] = 1 - 0.894119709
MSLBH2_cpt[:PGA=>:PGA_07, :MSLBH2=>:YES_MSLBH2] = 0.933362201
MSLBH2_cpt[:PGA=>:PGA_07, :MSLBH2=>:NO_MSLBH2] = 1 - 0.933362201
MSLBH2_cpt[:PGA=>:PGA_08, :MSLBH2=>:YES_MSLBH2] = 0.957292466
MSLBH2_cpt[:PGA=>:PGA_08, :MSLBH2=>:NO_MSLBH2] = 1 - 0.957292466
MSLBH2_cpt[:PGA=>:PGA_09, :MSLBH2=>:YES_MSLBH2] = 0.97213102
MSLBH2_cpt[:PGA=>:PGA_09, :MSLBH2=>:NO_MSLBH2] = 1 - 0.97213102
MSLBH2_cpt[:PGA=>:PGA_10, :MSLBH2=>:YES_MSLBH2] = 0.98149747
MSLBH2_cpt[:PGA=>:PGA_10, :MSLBH2=>:NO_MSLBH2] = 1 - 0.98149747
MSLBH2_cpt[:PGA=>:PGA_11, :MSLBH2=>:YES_MSLBH2] = 0.987515129
MSLBH2_cpt[:PGA=>:PGA_11, :MSLBH2=>:NO_MSLBH2] = 1 - 0.987515129
MSLBH2_cpt[:PGA=>:PGA_12, :MSLBH2=>:YES_MSLBH2] = 0.991447327
MSLBH2_cpt[:PGA=>:PGA_12, :MSLBH2=>:NO_MSLBH2] = 1 - 0.991447327
MSLBH2_cpt[:PGA=>:PGA_13, :MSLBH2=>:YES_MSLBH2] = 0.994058208
MSLBH2_cpt[:PGA=>:PGA_13, :MSLBH2=>:NO_MSLBH2] = 1 - 0.994058208
MSLBH2_cpt[:PGA=>:PGA_14, :MSLBH2=>:YES_MSLBH2] = 0.99581793
MSLBH2_cpt[:PGA=>:PGA_14, :MSLBH2=>:NO_MSLBH2] = 1 - 0.99581793
MSLBH2_cpt[:PGA=>:PGA_15, :MSLBH2=>:YES_MSLBH2] = 0.997020675
MSLBH2_cpt[:PGA=>:PGA_15, :MSLBH2=>:NO_MSLBH2] = 1 - 0.997020675
MSLBH2_cpt[:PGA=>:PGA_16, :MSLBH2=>:YES_MSLBH2] = 0.99785352
MSLBH2_cpt[:PGA=>:PGA_16, :MSLBH2=>:NO_MSLBH2] = 1 - 0.99785352
MSLBH2_cpt[:PGA=>:PGA_17, :MSLBH2=>:YES_MSLBH2] = 0.99843728
MSLBH2_cpt[:PGA=>:PGA_17, :MSLBH2=>:NO_MSLBH2] = 1 - 0.99843728
MSLBH2_cpt[:PGA=>:PGA_18, :MSLBH2=>:YES_MSLBH2] = 0.998851119
MSLBH2_cpt[:PGA=>:PGA_18, :MSLBH2=>:NO_MSLBH2] = 1 - 0.998851119
MSLBH2_cpt[:PGA=>:PGA_19, :MSLBH2=>:YES_MSLBH2] = 0.999147624
MSLBH2_cpt[:PGA=>:PGA_19, :MSLBH2=>:NO_MSLBH2] = 1 - 0.999147624

MSLBH2_node = DiscreteNode(:MSLBH2, MSLBH2_cpt)

``` MSLBH2-TIME node
```
t_MSLBH2_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:MSLBH2)
t_MSLBH2_cpt[:MSLBH2=>:YES_MSLBH2] = Uniform(1, 1200)
t_MSLBH2_cpt[:MSLBH2=>:NO_MSLBH2] = Normal(1200, 0)
t_MSLBH2_node = ContinuousNode(:t_MSLBH2, t_MSLBH2_cpt)


``` EXPLOSION node
```
explosion_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}([:PGA, :MSLB, :EXP])
explosion_cpt[:PGA=>:PGA_00, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0
explosion_cpt[:PGA=>:PGA_00, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1
explosion_cpt[:PGA=>:PGA_00, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.646
explosion_cpt[:PGA=>:PGA_00, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.646
explosion_cpt[:PGA=>:PGA_01, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 1
explosion_cpt[:PGA=>:PGA_01, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 0
explosion_cpt[:PGA=>:PGA_01, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.646
explosion_cpt[:PGA=>:PGA_01, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.646
explosion_cpt[:PGA=>:PGA_02, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 1
explosion_cpt[:PGA=>:PGA_02, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 0
explosion_cpt[:PGA=>:PGA_02, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.646
explosion_cpt[:PGA=>:PGA_02, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.646
explosion_cpt[:PGA=>:PGA_03, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 1
explosion_cpt[:PGA=>:PGA_03, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 0
explosion_cpt[:PGA=>:PGA_03, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.646
explosion_cpt[:PGA=>:PGA_03, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.646
explosion_cpt[:PGA=>:PGA_04, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 1e-5
explosion_cpt[:PGA=>:PGA_04, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 1e-5
explosion_cpt[:PGA=>:PGA_04, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64601
explosion_cpt[:PGA=>:PGA_04, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64601
explosion_cpt[:PGA=>:PGA_05, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 1e-5
explosion_cpt[:PGA=>:PGA_05, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 1e-5
explosion_cpt[:PGA=>:PGA_05, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64601
explosion_cpt[:PGA=>:PGA_05, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64601
explosion_cpt[:PGA=>:PGA_06, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 5e-5
explosion_cpt[:PGA=>:PGA_06, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 5e-5
explosion_cpt[:PGA=>:PGA_06, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64605
explosion_cpt[:PGA=>:PGA_06, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64605
explosion_cpt[:PGA=>:PGA_07, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00017
explosion_cpt[:PGA=>:PGA_07, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00017
explosion_cpt[:PGA=>:PGA_07, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64617
explosion_cpt[:PGA=>:PGA_07, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64617
explosion_cpt[:PGA=>:PGA_08, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.0002
explosion_cpt[:PGA=>:PGA_08, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.0002
explosion_cpt[:PGA=>:PGA_08, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.6462
explosion_cpt[:PGA=>:PGA_08, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.6462
explosion_cpt[:PGA=>:PGA_09, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00038
explosion_cpt[:PGA=>:PGA_09, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00038
explosion_cpt[:PGA=>:PGA_09, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64638
explosion_cpt[:PGA=>:PGA_09, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64638
explosion_cpt[:PGA=>:PGA_10, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00052
explosion_cpt[:PGA=>:PGA_10, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00052
explosion_cpt[:PGA=>:PGA_10, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64652
explosion_cpt[:PGA=>:PGA_10, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64652
explosion_cpt[:PGA=>:PGA_11, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00107
explosion_cpt[:PGA=>:PGA_11, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00107
explosion_cpt[:PGA=>:PGA_11, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64707
explosion_cpt[:PGA=>:PGA_11, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64707
explosion_cpt[:PGA=>:PGA_12, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.0015
explosion_cpt[:PGA=>:PGA_12, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.0015
explosion_cpt[:PGA=>:PGA_12, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.6475
explosion_cpt[:PGA=>:PGA_12, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.6475
explosion_cpt[:PGA=>:PGA_13, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00212
explosion_cpt[:PGA=>:PGA_13, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00212
explosion_cpt[:PGA=>:PGA_13, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64812
explosion_cpt[:PGA=>:PGA_13, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64812
explosion_cpt[:PGA=>:PGA_14, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00327
explosion_cpt[:PGA=>:PGA_14, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00327
explosion_cpt[:PGA=>:PGA_14, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.64927
explosion_cpt[:PGA=>:PGA_14, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.64927
explosion_cpt[:PGA=>:PGA_15, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00435
explosion_cpt[:PGA=>:PGA_15, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00435
explosion_cpt[:PGA=>:PGA_15, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.65035
explosion_cpt[:PGA=>:PGA_15, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.65035
explosion_cpt[:PGA=>:PGA_16, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00597
explosion_cpt[:PGA=>:PGA_16, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00597
explosion_cpt[:PGA=>:PGA_16, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.65197
explosion_cpt[:PGA=>:PGA_16, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.65197
explosion_cpt[:PGA=>:PGA_17, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.00798
explosion_cpt[:PGA=>:PGA_17, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.00798
explosion_cpt[:PGA=>:PGA_17, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.65398
explosion_cpt[:PGA=>:PGA_17, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.65398
explosion_cpt[:PGA=>:PGA_18, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.01082
explosion_cpt[:PGA=>:PGA_18, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.01082
explosion_cpt[:PGA=>:PGA_18, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.65682
explosion_cpt[:PGA=>:PGA_18, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.65682
explosion_cpt[:PGA=>:PGA_19, :MSLB=>:NO_MSLB, :EXP=>:YES_EXP] = 0.01584
explosion_cpt[:PGA=>:PGA_19, :MSLB=>:NO_MSLB, :EXP=>:NO_EXP] = 1 - 0.01584
explosion_cpt[:PGA=>:PGA_19, :MSLB=>:YES_MSLB, :EXP=>:YES_EXP] = 0.66184
explosion_cpt[:PGA=>:PGA_19, :MSLB=>:YES_MSLB, :EXP=>:NO_EXP] = 1 - 0.66184

explosion_node = DiscreteNode(:EXP, explosion_cpt)


``` LOOPH2 node
```
LOOPH2_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}([:EXP, :LOOPH2])
LOOPH2_cpt[:EXP=>:YES_EXP, :LOOPH2=>:YES_LOOPH2] = 0.004
LOOPH2_cpt[:EXP=>:YES_EXP, :LOOPH2=>:NO_LOOPH2] = 1 - 0.004
LOOPH2_cpt[:EXP=>:NO_EXP, :LOOPH2=>:YES_LOOPH2] = 0
LOOPH2_cpt[:EXP=>:NO_EXP, :LOOPH2=>:NO_LOOPH2] = 1
LOOPH2_node = DiscreteNode(:LOOPH2, LOOPH2_cpt)


``` LOOPH2-TIME node
```
t_LOOPH2_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LOOPH2)
t_LOOPH2_cpt[:LOOPH2=>:YES_LOOPH2] = Uniform(1, 1200)
t_LOOPH2_cpt[:LOOPH2=>:NO_LOOPH2] = Normal(1200, 0)
t_LOOPH2_node = ContinuousNode(:t_LOOPH2, t_LOOPH2_cpt)


``` OC node
```
OC_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}([:PGA, :OC])
OC_cpt[:PGA=>:PGA_00, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_00, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_01, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_01, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_02, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_02, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_03, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_03, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_04, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_04, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_05, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_05, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_06, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_06, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_07, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_07, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_08, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_08, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_09, :OC=>:YES_OC] = 0
OC_cpt[:PGA=>:PGA_09, :OC=>:NO_OC] = 1
OC_cpt[:PGA=>:PGA_10, :OC=>:YES_OC] = 4e-05
OC_cpt[:PGA=>:PGA_10, :OC=>:NO_OC] = 1 - 4e-05
OC_cpt[:PGA=>:PGA_11, :OC=>:YES_OC] = 0.00018
OC_cpt[:PGA=>:PGA_11, :OC=>:NO_OC] = 1 - 0.00018
OC_cpt[:PGA=>:PGA_12, :OC=>:YES_OC] = 0.00032
OC_cpt[:PGA=>:PGA_12, :OC=>:NO_OC] = 1 - 0.00032
OC_cpt[:PGA=>:PGA_13, :OC=>:YES_OC] = 0.00072
OC_cpt[:PGA=>:PGA_13, :OC=>:NO_OC] = 1 - 0.00072
OC_cpt[:PGA=>:PGA_14, :OC=>:YES_OC] = 0.00107
OC_cpt[:PGA=>:PGA_14, :OC=>:NO_OC] = 1 - 0.00107
OC_cpt[:PGA=>:PGA_15, :OC=>:YES_OC] = 0.00216
OC_cpt[:PGA=>:PGA_15, :OC=>:NO_OC] = 1 - 0.00216
OC_cpt[:PGA=>:PGA_16, :OC=>:YES_OC] = 0.00279
OC_cpt[:PGA=>:PGA_16, :OC=>:NO_OC] = 1 - 0.00279
OC_cpt[:PGA=>:PGA_17, :OC=>:YES_OC] = 0.00604
OC_cpt[:PGA=>:PGA_17, :OC=>:NO_OC] = 1 - 0.00604
OC_cpt[:PGA=>:PGA_18, :OC=>:YES_OC] = 0.01167
OC_cpt[:PGA=>:PGA_18, :OC=>:NO_OC] = 1 - 0.01167
OC_cpt[:PGA=>:PGA_19, :OC=>:YES_OC] = 0.01938
OC_cpt[:PGA=>:PGA_19, :OC=>:NO_OC] = 1 - 0.01938

OC_node = DiscreteNode(:OC, OC_cpt)


``` LHS node
```
LHS_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}([:OC, :EXP, :LHS])
LHS_cpt[:OC=>:YES_OC, :EXP=>:YES_EXP, :LHS=>:YES_LHS] = 1
LHS_cpt[:OC=>:YES_OC, :EXP=>:YES_EXP, :LHS=>:NO_LHS] = 0
LHS_cpt[:OC=>:YES_OC, :EXP=>:NO_EXP, :LHS=>:YES_LHS] = 0.0013
LHS_cpt[:OC=>:YES_OC, :EXP=>:NO_EXP, :LHS=>:NO_LHS] = 1 - 0.0013
LHS_cpt[:OC=>:NO_OC, :EXP=>:YES_EXP, :LHS=>:YES_LHS] = 1
LHS_cpt[:OC=>:NO_OC, :EXP=>:YES_EXP, :LHS=>:NO_LHS] = 0
LHS_cpt[:OC=>:NO_OC, :EXP=>:NO_EXP, :LHS=>:YES_LHS] = 0
LHS_cpt[:OC=>:NO_OC, :EXP=>:NO_EXP, :LHS=>:NO_LHS] = 1
LHS_node = DiscreteNode(:LHS, LHS_cpt)


``` LHS-TIME node
```
t_LHS_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:LHS)
t_LHS_cpt[:LHS=>:YES_LHS] = Uniform(1, 1200)
t_LHS_cpt[:LHS=>:NO_LHS] = Normal(1200, 0)
t_LHS_node = ContinuousNode(:t_LHS, t_LHS_cpt)


``` ACS node
```
data = CSV.read("networks/csv/ACS_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
rename!(data, :AGE_0 => :AGE_00)
df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
df_long2 = deepcopy(df_long1)
df_long2[!, :Π] = 1 .- df_long2[!, :Π]
df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
insertcols!(df_long1, 3, :ACS => fill(:YES_ACS, nrow(df_long1)))
insertcols!(df_long2, 3, :ACS => fill(:NO_ACS, nrow(df_long2)))
df_acs = sort(vcat(df_long1, df_long2))
acs_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_acs)
acs_node = DiscreteNode(:ACS, acs_cpt)


``` ACS-TIME node
```
t_acs_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:ACS)
t_acs_cpt[:ACS=>:YES_ACS] = Uniform(1, 1200)
t_acs_cpt[:ACS=>:NO_ACS] = Normal(1200, 0)
t_acs_node = ContinuousNode(:t_acs, t_acs_cpt)

``` ACS-rTIME node
```
rt_acs_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
rt_acs_cpt[] = Uniform(30.0, 90.0)
rt_acs_node = ContinuousNode(:rt_acs, rt_acs_cpt)


``` EDG node
```
data = CSV.read("networks/csv/EDG_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
rename!(data, :AGE_0 => :AGE_00)
df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
df_long2 = deepcopy(df_long1)
df_long2[!, :Π] = 1 .- df_long2[!, :Π]
df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
insertcols!(df_long1, 3, :EDG => fill(:YES_EDG, nrow(df_long1)))
insertcols!(df_long2, 3, :EDG => fill(:NO_EDG, nrow(df_long2)))
df_edg = sort(vcat(df_long1, df_long2))
edg_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_edg)
edg_node = DiscreteNode(:EDG, edg_cpt)


``` EDG-TIME node
```
t_edg_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:EDG)
t_edg_cpt[:EDG=>:YES_EDG] = Uniform(1, 1200)
t_edg_cpt[:EDG=>:NO_EDG] = Normal(1200, 0)
t_edg_node = ContinuousNode(:t_edg, t_edg_cpt)

``` EDG-rTIME node
```
rt_edg_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
rt_edg_cpt[] = Uniform(10.0, 180.0)
rt_edg_node = ContinuousNode(:rt_edg, rt_edg_cpt)


``` PDP node
```
data = CSV.read("networks/csv/PDP_probability.csv", DataFrame)
rename!(data, :Column1 => :PGA)
rename!(data, :AGE_0 => :AGE_00)
df_long1 = stack(data, Not(:PGA), variable_name=:AGE, value_name=:Π)
df_long2 = deepcopy(df_long1)
df_long2[!, :Π] = 1 .- df_long2[!, :Π]
df_long1[!, [:PGA, :AGE]] .= Symbol.(df_long1[!, [:PGA, :AGE]])
df_long2[!, [:PGA, :AGE]] .= Symbol.(df_long2[!, [:PGA, :AGE]])
insertcols!(df_long1, 3, :PDP => fill(:YES_PDP, nrow(df_long1)))
insertcols!(df_long2, 3, :PDP => fill(:NO_PDP, nrow(df_long2)))
df_pdp = sort(vcat(df_long1, df_long2))
pdp_cpt = DiscreteConditionalProbabilityTable{PreciseDiscreteProbability}(df_pdp)
pdp_node = DiscreteNode(:PDP, pdp_cpt)


``` PDP-TIME node
```
t_pdp_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}(:PDP)
t_pdp_cpt[:PDP=>:YES_PDP] = Uniform(1, 1200)
t_pdp_cpt[:PDP=>:NO_PDP] = Normal(1200, 0)
t_pdp_node = ContinuousNode(:t_pdp, t_pdp_cpt)

``` PDP-rTIME node
```
rt_pdp_cpt = ContinuousConditionalProbabilityTable{PreciseContinuousInput}()
rt_pdp_cpt[] = Uniform(10, 180)
rt_pdp_node = ContinuousNode(:rt_pdp, rt_pdp_cpt)


``` MODEL node
```
include("model_T.jl")

model_temp = ParallelModel(df -> model_temperatures(df.t_loca, df.t_loop, df.t_LHS, df.t_mslb, df.t_MSLBH2, df.t_LOOPH2, df.t_acs, df.rt_acs, df.t_edg, df.rt_edg, df.t_pdp, df.rt_pdp), :max_Ts)

function performance_function(threshold::Real, df::DataFrame)
    maxval = maximum(Matrix(df))
    return threshold - maxval
end

performance = df -> performance_function.(threshold, df.max_Ts)
sim = MonteCarlo(SIMULATIONS)
model_node = DiscreteFunctionalNode(:Reactor, [model_temp], performance, sim)


``` Enhanced Bayesian Networks
```
nodes = [
    age_node, pga_node,
    loca_node, t_loca_node,
    loop_node, t_loop_node,
    mslb_node, t_mslb_node,
    MSLBH2_node, t_MSLBH2_node,
    explosion_node,
    OC_node,
    LOOPH2_node, t_LOOPH2_node,
    LHS_node, t_LHS_node, acs_node, t_acs_node, rt_acs_node,
    edg_node, t_edg_node, rt_edg_node,
    pdp_node, t_pdp_node, rt_pdp_node,
    model_node
]

ebn = EnhancedBayesianNetwork(nodes)
add_child!(ebn, :AGE, :LOCA)
add_child!(ebn, :PGA, :LOCA)
add_child!(ebn, :LOCA, :t_loca)

add_child!(ebn, :PGA, :LOOP)
add_child!(ebn, :LOOP, :t_loop)

add_child!(ebn, :AGE, :MSLB)
add_child!(ebn, :PGA, :MSLB)
add_child!(ebn, :MSLB, :t_mslb)

add_child!(ebn, :PGA, :MSLBH2)
add_child!(ebn, :MSLBH2, :t_MSLBH2)
add_child!(ebn, :MSLB, :EXP)
add_child!(ebn, :PGA, :EXP)
add_child!(ebn, :EXP, :LOOPH2)
add_child!(ebn, :LOOPH2, :t_LOOPH2)

add_child!(ebn, :PGA, :OC)
add_child!(ebn, :EXP, :LHS)
add_child!(ebn, :OC, :LHS)
add_child!(ebn, :LHS, :t_LHS)

add_child!(ebn, :AGE, :ACS)
add_child!(ebn, :PGA, :ACS)
add_child!(ebn, :ACS, :t_acs)

add_child!(ebn, :AGE, :EDG)
add_child!(ebn, :PGA, :EDG)
add_child!(ebn, :EDG, :t_edg)

add_child!(ebn, :AGE, :PDP)
add_child!(ebn, :PGA, :PDP)
add_child!(ebn, :PDP, :t_pdp)

add_child!(ebn, :t_loca, :Reactor)
add_child!(ebn, :t_loop, :Reactor)
add_child!(ebn, :t_mslb, :Reactor)
add_child!(ebn, :t_MSLBH2, :Reactor)
add_child!(ebn, :t_LOOPH2, :Reactor)
add_child!(ebn, :t_LHS, :Reactor)
add_child!(ebn, :t_acs, :Reactor)
add_child!(ebn, :rt_acs, :Reactor)
add_child!(ebn, :t_edg, :Reactor)
add_child!(ebn, :rt_edg, :Reactor)
add_child!(ebn, :t_pdp, :Reactor)
add_child!(ebn, :rt_pdp, :Reactor)


order!(ebn)
gplot(ebn; NODESIZEFACTOR=0.1, ARROWLENGTH=0.05, NODELABELSIZE=2.5)
# --- initial time ---
t0 = time_ns()

# end

evaluate!(ebn, false, true)

# Save the network to disk
path_to_ebn = joinpath(current_dir, "networks", "ebn_jld2")
mkpath(path_to_ebn)
ebn_name = Dates.format(now(), "yyyy_mm_dd_HH_MM") * "_" *
           string(model_node.simulation) * ".jld2"
@save joinpath(path_to_ebn, ebn_name) ebn

# Print elapsed time
seconds = (time_ns() - t0) / 1e9
println("Elapsed time: $(round(seconds, digits=3)) s")

# rmprocs(workers())