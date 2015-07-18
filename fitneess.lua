Grid = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}

function calcFitness()
  local highNum = highestNum()
  local mono = monotonicity()
  local numFree = freeTiles()
end

function highestNum()
  local highest = -1
  for row in Grid do
    if math.max(row) > highest then
      highest = math.max(row)
    end
  end
  return math.log10(highest)
end

function monotonicity()

  local totals = {0,0,0,0}

  for x=1,5,1 do
    local current = 0
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
        nextValue = math.log10(Grid[x][current]) / math.log(2)
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

  for y=1,5,1 do
    local current = 0
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
end

function freeTiles():
  local count = 0
  for row in Grid do
    for val in row do
      if val == 0 then
        count = count + 1
      end
    end
  end
  return count
end

