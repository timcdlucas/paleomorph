riskr
======


An R package for calculating Risk battle odds.

![Rimmer from Red Dwarf](hqdefault.jpg)

> The only way was north; I had to go for it and pray the Gods were smiling on me. 
> I picked up the dice and threw two sixes. Caldecott couldn't believe it. 
> My go again; another two sixes!
>
> Arnold Rimmer


Rationale
----------

There are quite a few Risk battle odds calculators on the internet. 
But I didn't find one that calculates odds for a whole string of battles which is what we often want to know.
For example, you have 20 units in North Africa. 
What are the odds you can take Brazil (2 units), Venezuala (3 units) and Central America (10 units) and thus deny your opponent their N. America bonus reinforcements?

    1> fight(a = 10, d = c(2,3,10))
    Attacker wins 8%


Basic usage
------------

To install:

    library(devtools)
    install_github('timcdlucas/riskr')
    library(riskr)

Then the main function is `fight`. 
Give the number of attacking units (remember to minus one as you have to leave one behind). 
Also set the number of defensive units. This can be a single number, or a vector of numbers for a sequence of countries that you wish to take.

    fight(a = 10, d = 20)
    fight(a = 20, d = c(10, 1, 1, 1))
    fight(a = 20, d = 20, sims = 2000)


For more details, save the object and either plot it or get a short summary by printing.

    > f <- fight(a = 20, d = c(3, 4, 1, 10))
    Attacker wins 60.7%
    Attacker average loses: 16.5
    
    > f
    Risk battle odds for 20 attackers and
     18 defenders in 5 territories:

    Attacker survives 60.7%
    Attacker average loses: 16.5


    From 1000 simulations.
    > plot(f)

    



Next steps
------------

I realise this is up there with the pretty nerdy stuff I've done. But if there's rules (newer risk?) I've missed let me know. 
If you want odds for other games, also let me know and there's a vague change I might one day add it.

If you would like to add any thing, pull and push away.

To do: 
- Write a summary function that tells you how many units the attacker has left or how many countries they took if all their units died.
- Improve histograms.


