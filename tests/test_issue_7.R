library(climdex.pcic)
library(RUnit)

climdex.pcic.test.one.column.dates <- function() {
    ci <- climdexInput.csv(tmax.file='single_date_field_tmax.csv',
                           tmin.file='single_date_field_tmin.csv',
                           prec.file='single_date_field_prec.csv',
                           date.types=list(list(fields=c("date"), format="%Y-%j")))
    checkTrue(inherits(ci, 'climdexInput'))
}
