zpoly <- ZonalPol(n = 4, lambda = c(2, 2))
compactSymmetricQspray(zpoly)

schurCombo <- SchurCombination(zpoly)

result <- qzero()
for(lst2 in schurCombo) {
  result <- result + lst2[["coeff"]] * SchurPol(n = 4, lambda = lst2[["lambda"]])
}

print(result == zpoly)
