def presentffuture(i, n, FV=1):
    ''' Future value to Present Value
    where FV = future value of money, i is the interest rate etc, n is the number of years. '''
    return FV * (1 + i) ** -n


def futurefpresent(i, n, PV=1):
    ''' Present Value to Future Value
    where FV = future value of money, i is the interest rate etc, n is the number of years. '''
    return PV * (1 + i) ** n


def futurefannual(i, n, A=1):
    return A * ((1 + i) ** n - 1) / i


def annual2present(i, n, A=1):
    return A * ((1 + i) ** n - 1) / (i * (1 + i) ** n)


def annual2present_geometric(i, n, g, A1=1):
    return A1 * (1 - (1 + g) ** n * (1 + i) ** -n) / (i - g)


def present2annual(i, n, PV=1):
    '''Present Value to Annualised'''
    return PV * (i * (1 + i) ** n) / ((1 + i) ** n - 1)