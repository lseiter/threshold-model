library(RSQLite)
db <- dbConnect(SQLite(), dbname="voting.sqlite")

dbSendQuery(db, "DROP TABLE if exists session")
dbSendQuery(db, "CREATE TABLE session(
           sessionid INTEGER PRIMARY KEY, -- Autoincrement
           configuration text,
          iterations integer,
           k integer,
           type text,
           distanceNeighbors integer,
           probRewire real,  
           n integer,
           step integer,
           result text,
           finalPercentages text,
           seed text,
           tb text,
           decay real,
           winner integer,
           winnerPercent real,
           winnerLead  real,
           meanDistance real, 
           clusterCoeff real,
           meanDegree real,
           degreeDistribution text
)")

