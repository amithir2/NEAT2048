-- MarI/O by SethBling
-- Feel free to use this code, but please do not redistribute it.
-- Intended for use with the BizHawk emulator and Super Mario World or Super Mario Bros. ROM.
-- For SMW, make sure you have a save state named "DP1.state" at the beginning of a level,
-- and put a copy in both the Lua folder and the root directory of BizHawk.
 
-- removed variable Filename
-- what is Filename?
ButtonNames = {
        "Up",
        "Down",
        "Left",
        "Right",
}
 
BoxRadius = 1
Inputs = 16
Outputs = #ButtonNames

Grid = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}

Population = 300
DeltaDisjoint = 2.0
DeltaWeights = 0.4
DeltaThreshold = 1.0
 
StaleSpecies = 15
 
MutateConnectionsChance = 0.25
PerturbChance = 0.90
CrossoverChance = 0.75
LinkMutationChance = 2.0
NodeMutationChance = 0.50
BiasMutationChance = 0.40
StepSize = 0.1
DisableMutationChance = 0.4
EnableMutationChance = 0.2
 
MaxNodes = 1000000

function printArray()
        print(table.concat(Grid[1], "\t"))
        print(table.concat(Grid[2], "\t"))
        print(table.concat(Grid[3], "\t"))
        print(table.concat(Grid[4], "\t"))
end

function simulateClear()
        Grid = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}
end

function max(t, fn)
    if #t == 0 then return nil end
    local key, value = 1, t[1]
    for i = 2, #t do
        if fn(value, t[i]) then
            key, value = i, t[i]
        end
    end
    return value
end

function simulateUpdate() 
        for y=1,4,1 do 
                for x=1,4,1 do
                        if Grid[y][x] == 0 then
                                Grid[y][x] = 2
                                return true
                        end
                end
        end
        return false

end

function moveToLeft(row)
        for i=4,1,-1 do
                if row[i] == 0 then 
                        table.remove(row, i)
                end
        end
        for i=1,#row-1,1 do 
                local j = i+1
                if row[i] == row[j] then
                        row[i] = row[i]*2
                        row[j] = 0
                end
        end
        for i=4,1,-1 do
                if row[i] == 0 then 
                        table.remove(row, i)
                end
        end
        while #row ~= 4 do
                table.insert(row, #row+1, 0)
        end
        return row
end

function moveToRight(row)
        for i=4,1,-1 do
                if row[i] == 0 then 
                        table.remove(row, i)
                end
        end
        for i=#row,2,-1 do 
                local j = i-1
                if row[i] == row[j] then
                        row[i] = row[i]*2
                        row[j] = 0
                end
        end
        for i=4,1,-1 do
                if row[i] == 0 then 
                        table.remove(row, i)
                end
        end
        while #row ~= 4 do
                table.insert(row, 1, 0)
        end
        return row
end

function simulateSendResults(controller) 
        if controller["P1 Left"] then
                Grid[1] = moveToLeft(Grid[1])
                Grid[2] = moveToLeft(Grid[2])
                Grid[3] = moveToLeft(Grid[3])
                Grid[4] = moveToLeft(Grid[4])
                return
        end
        if controller["P1 Right"] then
                Grid[1] = moveToRight(Grid[1])
                Grid[2] = moveToRight(Grid[2])
                Grid[3] = moveToRight(Grid[3])
                Grid[4] = moveToRight(Grid[4])
                return
        end
        if controller["P1 Up"] then
                -- print("start")
                -- printArray()
                col1 = {Grid[1][1], Grid[2][1], Grid[3][1], Grid[4][1]}
                col2 = {Grid[1][2], Grid[2][2], Grid[3][2], Grid[4][2]}
                col3 = {Grid[1][3], Grid[2][3], Grid[3][3], Grid[4][3]}
                col4 = {Grid[1][4], Grid[2][4], Grid[3][4], Grid[4][4]}
                col1 = moveToLeft(col1)
                col2 = moveToLeft(col2)
                col3 = moveToLeft(col3)
                col4 = moveToLeft(col4)
                
                simulateClear()
                for i=1,4,1 do
                        Grid[i][1] = col1[i]
                        Grid[i][2] = col2[i]
                        Grid[i][3] = col3[i]
                        Grid[i][4] = col4[i]
                end
                -- printArray()
                -- print("stop")
                
                return
        end
        if controller["P1 Down"] then
                -- print("start")
                -- printArray()
                col1 = {Grid[1][1], Grid[2][1], Grid[3][1], Grid[4][1]}
                col2 = {Grid[1][2], Grid[2][2], Grid[3][2], Grid[4][2]}
                col3 = {Grid[1][3], Grid[2][3], Grid[3][3], Grid[4][3]}
                col4 = {Grid[1][4], Grid[2][4], Grid[3][4], Grid[4][4]}
                col1 = moveToRight(col1)
                col2 = moveToRight(col2)
                col3 = moveToRight(col3)
                col4 = moveToRight(col4)
                simulateClear()
                for i=1,4,1 do
                        Grid[i][1] = col1[i]
                        Grid[i][2] = col2[i]
                        Grid[i][3] = col3[i]
                        Grid[i][4] = col4[i]
                end
                -- printArray()
                -- print("stop")
                return
        end

end

function getTile(y, x)
        if math.min(x, y) < 0 or math.max(x, y) >= 4 then 
                return 0
        end
        return Grid[y][x]
end
 
function getInputs()
       
        local inputs = {}
       
        for y=1,4,1 do
                for x=1,4,1 do
                        inputs[#inputs+1] = getTile(y, x)
                end
        end 
       
        return inputs
end
 
function sigmoid(x)
        return 2/(1+math.exp(-4.9*x))-1
end
 
function newInnovation()
        pool.innovation = pool.innovation + 1
        return pool.innovation
end
 
function newPool()
        local pool = {}
        pool.species = {}
        pool.generation = 0
        pool.innovation = Outputs
        pool.currentSpecies = 1
        pool.currentGenome = 1
        pool.currentFrame = 0
        pool.maxFitness = 0
       
        return pool
end
 
function newSpecies()
        local species = {}
        species.topFitness = 0
        species.staleness = 0
        species.genomes = {}
        species.averageFitness = 0
       
        return species
end
 
function newGenome()
        local genome = {}
        genome.genes = {}
        genome.fitness = 0
        genome.adjustedFitness = 0
        genome.network = {}
        genome.maxneuron = 0
        genome.globalRank = 0
        genome.mutationRates = {}
        genome.mutationRates["connections"] = MutateConnectionsChance
        genome.mutationRates["link"] = LinkMutationChance
        genome.mutationRates["bias"] = BiasMutationChance
        genome.mutationRates["node"] = NodeMutationChance
        genome.mutationRates["enable"] = EnableMutationChance
        genome.mutationRates["disable"] = DisableMutationChance
        genome.mutationRates["step"] = StepSize
       
        return genome
end
 
function copyGenome(genome)
        local genome2 = newGenome()
        for g=1,#genome.genes do
                table.insert(genome2.genes, copyGene(genome.genes[g]))
        end
        genome2.maxneuron = genome.maxneuron
        genome2.mutationRates["connections"] = genome.mutationRates["connections"]
        genome2.mutationRates["link"] = genome.mutationRates["link"]
        genome2.mutationRates["bias"] = genome.mutationRates["bias"]
        genome2.mutationRates["node"] = genome.mutationRates["node"]
        genome2.mutationRates["enable"] = genome.mutationRates["enable"]
        genome2.mutationRates["disable"] = genome.mutationRates["disable"]
       
        return genome2
end
 
function basicGenome()
        local genome = newGenome()
        local innovation = 1
 
        genome.maxneuron = Inputs
        mutate(genome)
       
        return genome
end
 
function newGene()
        local gene = {}
        gene.into = 0
        gene.out = 0
        gene.weight = 0.0
        gene.enabled = true
        gene.innovation = 0
       
        return gene
end
 
function copyGene(gene)
        local gene2 = newGene()
        gene2.into = gene.into
        gene2.out = gene.out
        gene2.weight = gene.weight
        gene2.enabled = gene.enabled
        gene2.innovation = gene.innovation
       
        return gene2
end
 
function newNeuron()
        local neuron = {}
        neuron.incoming = {}
        neuron.value = 0.0
       
        return neuron
end
 
function generateNetwork(genome)
        local network = {}
        network.neurons = {}
       
        for i=1,Inputs do
                network.neurons[i] = newNeuron()
        end
       
        for o=1,Outputs do
                network.neurons[MaxNodes+o] = newNeuron()
        end
       
        table.sort(genome.genes, function (a,b)
                return (a.out < b.out)
        end)
        for i=1,#genome.genes do
                local gene = genome.genes[i]
                if gene.enabled then
                        if network.neurons[gene.out] == nil then
                                network.neurons[gene.out] = newNeuron()
                        end
                        local neuron = network.neurons[gene.out]
                        table.insert(neuron.incoming, gene)
                        if network.neurons[gene.into] == nil then
                                network.neurons[gene.into] = newNeuron()
                        end
                end
        end
       
        genome.network = network
end
 
 -- update this logic to 2048 movement logic
function evaluateNetwork(network, inputs)
        if #inputs ~= Inputs then
                print("Incorrect number of neural network inputs.")
                return {}
        end
       
        for i=1,Inputs do
                network.neurons[i].value = inputs[i]
        end
       
        for _,neuron in pairs(network.neurons) do
                local sum = 0
                for j = 1,#neuron.incoming do
                        local incoming = neuron.incoming[j]
                        local other = network.neurons[incoming.into]
                        sum = sum + incoming.weight * other.value
                end
               
                if #neuron.incoming > 0 then
                        neuron.value = sigmoid(sum)
                end
        end
       
        local outputs = {}
        for o=1,Outputs do
                local button = "P1 " .. ButtonNames[o]
                if network.neurons[MaxNodes+o].value > 0 then
                        outputs[button] = true
                else
                        outputs[button] = false
                end
        end
       
        return outputs
end
 
function crossover(g1, g2)
        -- Make sure g1 is the higher fitness genome
        if g2.fitness > g1.fitness then
                tempg = g1
                g1 = g2
                g2 = tempg
        end
 
        local child = newGenome()
       
        local innovations2 = {}
        for i=1,#g2.genes do
                local gene = g2.genes[i]
                innovations2[gene.innovation] = gene
        end
       
        for i=1,#g1.genes do
                local gene1 = g1.genes[i]
                local gene2 = innovations2[gene1.innovation]
                if gene2 ~= nil and math.random(2) == 1 and gene2.enabled then
                        table.insert(child.genes, copyGene(gene2))
                else
                        table.insert(child.genes, copyGene(gene1))
                end
        end
       
        child.maxneuron = math.max(g1.maxneuron,g2.maxneuron)
       
        for mutation,rate in pairs(g1.mutationRates) do
                child.mutationRates[mutation] = rate
        end
       
        return child
end
 
function randomNeuron(genes, nonInput)
        local neurons = {}
        if not nonInput then
                for i=1,Inputs do
                        neurons[i] = true
                end
        end
        for o=1,Outputs do
                neurons[MaxNodes+o] = true
        end
        for i=1,#genes do
                if (not nonInput) or genes[i].into > Inputs then
                        neurons[genes[i].into] = true
                end
                if (not nonInput) or genes[i].out > Inputs then
                        neurons[genes[i].out] = true
                end
        end
 
        local count = 0
        for _,_ in pairs(neurons) do
                count = count + 1
        end
        local n = math.random(1, count)
       
        for k,v in pairs(neurons) do
                n = n-1
                if n == 0 then
                        return k
                end
        end
       
        return 0
end
 
function containsLink(genes, link)
        for i=1,#genes do
                local gene = genes[i]
                if gene.into == link.into and gene.out == link.out then
                        return true
                end
        end
end
 
function pointMutate(genome)
        local step = genome.mutationRates["step"]
       
        for i=1,#genome.genes do
                local gene = genome.genes[i]
                if math.random() < PerturbChance then
                        gene.weight = gene.weight + math.random() * step*2 - step
                else
                        gene.weight = math.random()*4-2
                end
        end
end
 
function linkMutate(genome, forceBias)
        local neuron1 = randomNeuron(genome.genes, false)
        local neuron2 = randomNeuron(genome.genes, true)
         
        local newLink = newGene()
        if neuron1 <= Inputs and neuron2 <= Inputs then
                --Both input nodes
                return
        end
        if neuron2 <= Inputs then
                -- Swap output and input
                local temp = neuron1
                neuron1 = neuron2
                neuron2 = temp
        end
 
        newLink.into = neuron1
        newLink.out = neuron2
        if forceBias then
                newLink.into = Inputs
        end
       
        if containsLink(genome.genes, newLink) then
                return
        end
        newLink.innovation = newInnovation()
        newLink.weight = math.random()*4-2
       
        table.insert(genome.genes, newLink)
end
 
function nodeMutate(genome)
        if #genome.genes == 0 then
                return
        end
 
        genome.maxneuron = genome.maxneuron + 1
 
        local gene = genome.genes[math.random(1,#genome.genes)]
        if not gene.enabled then
                return
        end
        gene.enabled = false
       
        local gene1 = copyGene(gene)
        gene1.out = genome.maxneuron
        gene1.weight = 1.0
        gene1.innovation = newInnovation()
        gene1.enabled = true
        table.insert(genome.genes, gene1)
       
        local gene2 = copyGene(gene)
        gene2.into = genome.maxneuron
        gene2.innovation = newInnovation()
        gene2.enabled = true
        table.insert(genome.genes, gene2)
end
 
function enableDisableMutate(genome, enable)
        local candidates = {}
        for _,gene in pairs(genome.genes) do
                if gene.enabled == not enable then
                        table.insert(candidates, gene)
                end
        end
       
        if #candidates == 0 then
                return
        end
       
        local gene = candidates[math.random(1,#candidates)]
        gene.enabled = not gene.enabled
end
 
function mutate(genome)
        for mutation,rate in pairs(genome.mutationRates) do
                if math.random(1,2) == 1 then
                        genome.mutationRates[mutation] = 0.95*rate
                else
                        genome.mutationRates[mutation] = 1.05263*rate
                end
        end
 
        if math.random() < genome.mutationRates["connections"] then
                pointMutate(genome)
        end
       
        local p = genome.mutationRates["link"]
        while p > 0 do
                if math.random() < p then
                        linkMutate(genome, false)
                end
                p = p - 1
        end
 
        p = genome.mutationRates["bias"]
        while p > 0 do
                if math.random() < p then
                        linkMutate(genome, true)
                end
                p = p - 1
        end
       
        p = genome.mutationRates["node"]
        while p > 0 do
                if math.random() < p then
                        nodeMutate(genome)
                end
                p = p - 1
        end
       
        p = genome.mutationRates["enable"]
        while p > 0 do
                if math.random() < p then
                        enableDisableMutate(genome, true)
                end
                p = p - 1
        end
 
        p = genome.mutationRates["disable"]
        while p > 0 do
                if math.random() < p then
                        enableDisableMutate(genome, false)
                end
                p = p - 1
        end
end
 
function disjoint(genes1, genes2)
        local i1 = {}
        for i = 1,#genes1 do
                local gene = genes1[i]
                i1[gene.innovation] = true
        end
 
        local i2 = {}
        for i = 1,#genes2 do
                local gene = genes2[i]
                i2[gene.innovation] = true
        end
       
        local disjointGenes = 0
        for i = 1,#genes1 do
                local gene = genes1[i]
                if not i2[gene.innovation] then
                        disjointGenes = disjointGenes+1
                end
        end
       
        for i = 1,#genes2 do
                local gene = genes2[i]
                if not i1[gene.innovation] then
                        disjointGenes = disjointGenes+1
                end
        end
       
        local n = math.max(#genes1, #genes2)
       
        return disjointGenes / n
end
 
function weights(genes1, genes2)
        local i2 = {}
        for i = 1,#genes2 do
                local gene = genes2[i]
                i2[gene.innovation] = gene
        end
 
        local sum = 0
        local coincident = 0
        for i = 1,#genes1 do
                local gene = genes1[i]
                if i2[gene.innovation] ~= nil then
                        local gene2 = i2[gene.innovation]
                        sum = sum + math.abs(gene.weight - gene2.weight)
                        coincident = coincident + 1
                end
        end
       
        return sum / coincident
end
       
function sameSpecies(genome1, genome2)
        local dd = DeltaDisjoint*disjoint(genome1.genes, genome2.genes)
        local dw = DeltaWeights*weights(genome1.genes, genome2.genes)
        return dd + dw < DeltaThreshold
end
 
function rankGlobally()
        local global = {}
        for s = 1,#pool.species do
                local species = pool.species[s]
                for g = 1,#species.genomes do
                        table.insert(global, species.genomes[g])
                end
        end
        table.sort(global, function (a,b)
                return (a.fitness < b.fitness)
        end)
       
        for g=1,#global do
                global[g].globalRank = g
        end
end
 
function calculateAverageFitness(species)
        local total = 0
       
        for g=1,#species.genomes do
                local genome = species.genomes[g]
                total = total + genome.globalRank
        end
       
        species.averageFitness = total / #species.genomes
end
 
function totalAverageFitness()
        local total = 0
        for s = 1,#pool.species do
                local species = pool.species[s]
                total = total + species.averageFitness
        end
 
        return total
end
 
function cullSpecies(cutToOne)
        for s = 1,#pool.species do
                local species = pool.species[s]
               
                table.sort(species.genomes, function (a,b)
                        return (a.fitness > b.fitness)
                end)
               
                local remaining = math.ceil(#species.genomes/2)
                if cutToOne then
                        remaining = 1
                end
                while #species.genomes > remaining do
                        table.remove(species.genomes)
                end
        end
end
 
function breedChild(species)
        local child = {}
        if math.random() < CrossoverChance then
                g1 = species.genomes[math.random(1, #species.genomes)]
                g2 = species.genomes[math.random(1, #species.genomes)]
                child = crossover(g1, g2)
        else
                g = species.genomes[math.random(1, #species.genomes)]
                child = copyGenome(g)
        end
       
        mutate(child)
       
        return child
end
 
function removeStaleSpecies()
        local survived = {}
 
        for s = 1,#pool.species do
                local species = pool.species[s]
               
                table.sort(species.genomes, function (a,b)
                        return (a.fitness > b.fitness)
                end)
               
                if species.genomes[1].fitness > species.topFitness then
                        species.topFitness = species.genomes[1].fitness
                        species.staleness = 0
                else
                        species.staleness = species.staleness + 1
                end
                if species.staleness < StaleSpecies or species.topFitness >= pool.maxFitness then
                        table.insert(survived, species)
                end
        end
 
        pool.species = survived
end
 
function removeWeakSpecies()
        local survived = {}
 
        local sum = totalAverageFitness()
        for s = 1,#pool.species do
                local species = pool.species[s]
                breed = math.floor(species.averageFitness / sum * Population)
                if breed >= 1 then
                        table.insert(survived, species)
                end
        end
 
        pool.species = survived
end
 
 
function addToSpecies(child)
        local foundSpecies = false
        for s=1,#pool.species do
                local species = pool.species[s]
                if not foundSpecies and sameSpecies(child, species.genomes[1]) then
                        table.insert(species.genomes, child)
                        foundSpecies = true
                end
        end
       
        if not foundSpecies then
                local childSpecies = newSpecies()
                table.insert(childSpecies.genomes, child)
                table.insert(pool.species, childSpecies)
        end
end
 
function newGeneration()
        cullSpecies(false) -- Cull the bottom half of each species
        rankGlobally()
        removeStaleSpecies()
        rankGlobally()
        for s = 1,#pool.species do
                local species = pool.species[s]
                calculateAverageFitness(species)
        end
        removeWeakSpecies()
        local sum = totalAverageFitness()
        local children = {}
        for s = 1,#pool.species do
                local species = pool.species[s]
                breed = math.floor(species.averageFitness / sum * Population) - 1
                for i=1,breed do
                        table.insert(children, breedChild(species))
                end
        end
        cullSpecies(true) -- Cull all but the top member of each species
        while #children + #pool.species < Population do
                local species = pool.species[math.random(1, #pool.species)]
                table.insert(children, breedChild(species))
        end
        for c=1,#children do
                local child = children[c]
                addToSpecies(child)
        end
       
        pool.generation = pool.generation + 1
       
        -- writeFile("backup." .. pool.generation .. "." .. forms.gettext(saveLoadFile))
end
       
function initializePool()
        pool = newPool()
 
        for i=1,Population do
                basic = basicGenome()
                addToSpecies(basic)
        end
 
        initializeRun()
end
 
function initializeRun()
        highestSquare = 0
        pool.currentFrame = 0

        -- clear results (don't send results to game)
        -- simulateClear()
       
        local species = pool.species[pool.currentSpecies]
        local genome = species.genomes[pool.currentGenome]

        -- generate network and neurons and etc
        generateNetwork(genome)
        -- evaluate valid current moves
        evaluateCurrent()
end
 
-- what is this logic even? 
function evaluateCurrent()
        local species = pool.species[pool.currentSpecies]
        local genome = species.genomes[pool.currentGenome]
 
        inputs = getInputs()
        controller = evaluateNetwork(genome.network, inputs)
        -- print(controll)
       
        -- if controller["P1 Left"] and controller["P1 Right"] then
        --         controller["P1 Left"] = false
        --         controller["P1 Right"] = false
        -- end
        -- if controller["P1 Up"] and controller["P1 Down"] then
        --         controller["P1 Up"] = false
        --         controller["P1 Down"] = false
        -- end
 
end
 
if pool == nil then
        initializePool()
end
 
 
function nextGenome()
        pool.currentGenome = pool.currentGenome + 1
        if pool.currentGenome > #pool.species[pool.currentSpecies].genomes then
                pool.currentGenome = 1
                pool.currentSpecies = pool.currentSpecies+1
                if pool.currentSpecies > #pool.species then
                        newGeneration()
                        pool.currentSpecies = 1
                end
        end
end
 
function fitnessAlreadyMeasured()
        local species = pool.species[pool.currentSpecies]
        local genome = species.genomes[pool.currentGenome]
       
        return genome.fitness ~= 0
end
 
function writeFile(filename)
        local file = io.open(filename, "w")
        file:write(pool.generation .. "\n")
        file:write(pool.maxFitness .. "\n")
        file:write(#pool.species .. "\n")
        for n,species in pairs(pool.species) do
                file:write(species.topFitness .. "\n")
                file:write(species.staleness .. "\n")
                file:write(#species.genomes .. "\n")
                for m,genome in pairs(species.genomes) do
                        file:write(genome.fitness .. "\n")
                        file:write(genome.maxneuron .. "\n")
                        for mutation,rate in pairs(genome.mutationRates) do
                                file:write(mutation .. "\n")
                                file:write(rate .. "\n")
                        end
                        file:write("done\n")
                       
                        file:write(#genome.genes .. "\n")
                        for l,gene in pairs(genome.genes) do
                                file:write(gene.into .. " ")
                                file:write(gene.out .. " ")
                                file:write(gene.weight .. " ")
                                file:write(gene.innovation .. " ")
                                if(gene.enabled) then
                                        file:write("1\n")
                                else
                                        file:write("0\n")
                                end
                        end
                end
        end
        file:close()
end
 
function savePool()
        local filename = forms.gettext(saveLoadFile)
        writeFile(filename)
end
 
function loadFile(filename)
        local file = io.open(filename, "r")
        pool = newPool()
        pool.generation = file:read("*number")
        pool.maxFitness = file:read("*number")
        forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
        local numSpecies = file:read("*number")
        for s=1,numSpecies do
                local species = newSpecies()
                table.insert(pool.species, species)
                species.topFitness = file:read("*number")
                species.staleness = file:read("*number")
                local numGenomes = file:read("*number")
                for g=1,numGenomes do
                        local genome = newGenome()
                        table.insert(species.genomes, genome)
                        genome.fitness = file:read("*number")
                        genome.maxneuron = file:read("*number")
                        local line = file:read("*line")
                        while line ~= "done" do
                                genome.mutationRates[line] = file:read("*number")
                                line = file:read("*line")
                        end
                        local numGenes = file:read("*number")
                        for n=1,numGenes do
                                local gene = newGene()
                                table.insert(genome.genes, gene)
                                local enabled
                                gene.into, gene.out, gene.weight, gene.innovation, enabled = file:read("*number", "*number", "*number", "*number", "*number")
                                if enabled == 0 then
                                        gene.enabled = false
                                else
                                        gene.enabled = true
                                end
                               
                        end
                end
        end
        file:close()
       
        while fitnessAlreadyMeasured() do
                nextGenome()
        end
        initializeRun()
        pool.currentFrame = pool.currentFrame + 1
end
 
function loadPool()
        local filename = forms.gettext(saveLoadFile)
        loadFile(filename)
end
 
function playTop()
        local maxfitness = 0
        local maxs, maxg
        for s,species in pairs(pool.species) do
                for g,genome in pairs(species.genomes) do
                        if genome.fitness > maxfitness then
                                maxfitness = genome.fitness
                                maxs = s
                                maxg = g
                        end
                end
        end
       
        pool.currentSpecies = maxs
        pool.currentGenome = maxg
        pool.maxFitness = maxfitness
        forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
        initializeRun()
        pool.currentFrame = pool.currentFrame + 1
        return
end
 
function onExit()
        forms.destroy(form)
end
 

 function calcFitness()
  local highNum = highestNum()
  local mono = monotonicity()
  local numFree = freeTiles()
  return highNum * 10 + mono + numFree

end

function highestNum()
  local highest = -1
  for i = 1,tablelength(Grid) do
    if max(Grid[i], function(a,b) return a < b end) > highest then
      highest = max(Grid[i], function(a,b) return a < b end)
    end
  end
  return math.log10(highest)
end

function monotonicity()

  local totals = {0,0,0,0}

  for x=1,4,1 do
    local current = 1
    local nextVal = current+1
    while nextVal < 4 do
      while nextVal < 4 and Grid[x][nextVal] == 0 do
        nextVal = nextVal + 1
      end
      if nextVal >= 4 then
        nextVal = nextVal - 1
      end
      local currentValue = -1
      if Grid[x][current] ~= 0 then
        currentValue = math.log10(Grid[x][current]) / math.log(2)
      end
      if Grid[x][current] == 0 then
        currentValue = 0
      end
      local nextValue = -1
      if Grid[x][nextVal] ~= 0 then
        nextValue = math.log10(Grid[x][nextVal]) / math.log(2)
      end
      if Grid[x][nextVal] == 0 then
        nextValue = 0
      end
      if currentValue > nextValue then
        totals[1] = totals[1] + nextValue - currentValue
      end
      if nextValue > currentValue then
        totals[2] = totals[2] + currentValue - nextValue
      end
      current = nextVal
      nextVal = nextVal + 1
    end
  end

  for y=1,4,1 do
    local current = 1
    local nextVal = current+1
    while nextVal < 4 do
      while nextVal < 4 and Grid[nextVal][y] == 0 do
        nextVal = nextVal + 1
      end
      if nextVal >= 4 then
        nextVal = nextVal - 1
      end
      local currentValue = -1
      if Grid[current][y] ~= 0 then
        currentValue = math.log10(Grid[current][y]) / math.log(2)
      end
      if Grid[current][y] == 0 then
        currentValue = 0
      end
      local nextValue = -1
      if Grid[current][y] ~= 0 then
        nextValue = math.log10(Grid[current][y]) / math.log(2)
      end
      if Grid[current][y] == 0 then
        nextValue = 0
      end
      if currentValue > nextValue then
        totals[3] = totals[3] + nextValue - currentValue
      end
      if nextValue > currentValue then
        totals[4] = totals[4] + currentValue - nextValue
      end
      current = nextVal
      nextVal = nextVal + 1
    end
  end

  return math.max(totals[1], totals[2]) + math.max(totals[3], totals[4])

end

function freeTiles()
  local count = 0
  for i = 1,tablelength(Grid) do
    for j=1,tablelength(Grid[i]) do
      if j == 0 then
        count = count + 1
      end
    end
  end
  return count
end

function tablelength(T)
  local count = 0
  for _ in pairs(T) do count = count + 1 end
  return count
end


writeFile("temp.pool")
 
 -- TODO add an exit event
-- event.onexit(onExit)
 

-- prevGrid = Grid
stillPlaying = true
bestScore = 0
numberOfMoves = 0
while stillPlaying do
        numberOfMoves = numberOfMoves + 1

        local species = pool.species[pool.currentSpecies]
        local genome = species.genomes[pool.currentGenome]

        -- this evaluates the current frame and determines possible outputs
        evaluateCurrent()
        
        -- need to send results to game
        simulateSendResults(controller)
        -- printArray()



        row1 = max(Grid[1], function(a,b) return a < b end)
        row2 = max(Grid[2], function(a,b) return a < b end)
        row3 = max(Grid[3], function(a,b) return a < b end)
        row4 = max(Grid[4], function(a,b) return a < b end)
        bestScore = math.max(row1, row2, row3, row4) 


        -- our fitness heuristic
        local fitness = calcFitness()
        -- print("Max Fitness: " .. fitness)
        if fitness == 0 then
                fitness = -1
        end
        genome.fitness = fitness



       
        if fitness > pool.maxFitness then
                pool.maxFitness = fitness
                -- 
        end
       
       -- print("Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " fitness: " .. fitness)
        pool.currentSpecies = 1
        pool.currentGenome = 1

        -- while fitness not equal to 0, keep iterating
        while fitnessAlreadyMeasured() do
                nextGenome()
        end
        initializeRun()
 
        local measured = 0
        local total = 0
        for _,species in pairs(pool.species) do
                for _,genome in pairs(species.genomes) do
                        total = total + 1
                        if genome.fitness ~= 0 then
                                measured = measured + 1
                        end
                end
        end
        -- print("Gen " .. pool.generation .. " species " .. pool.currentSpecies .. " genome " .. pool.currentGenome .. " (" .. math.floor(measured/total*100) .. "%)")
        -- print("Fitness: " .. math.floor(highestSquare - (pool.currentFrame) / 2))
        -- print("Max Fitness: " .. math.floor(pool.maxFitness))
               
        pool.currentFrame = pool.currentFrame + 1
 
        -- iterate to next frame somehow
        stillPlaying = simulateUpdate()
        if stillPlaying == false then
                print("Score: " .. bestScore .. " Number of Moves: " .. numberOfMoves)
                simulateClear()
                numberOfMoves = 0
                stillPlaying = true
        end
end