Grid = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}

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

function max(t, fn)
    if #t == 0 then return nil, nil end
    local key, value = 1, t[1]
    for i = 2, #t do
        if fn(value, t[i]) then
            key, value = i, t[i]
        end
    end
    return key, value
end

calcFitness()
