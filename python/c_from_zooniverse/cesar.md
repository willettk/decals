

Cesar has two outputs.

'Extracts' records each classification by each user, just like the script get_subject_classifications.
'Reductions' aggregates these classifications (for now, by adding them up).

In both cases, the `data.X` columns are the answers' index value (e.g. 0=frist option, 1=second option, etc...).  
For 'reductions', the value is how many users voted for that answer.
For 'extracts', the value is presumably always 1.

You can access Cesar data with urls like this:
 [https://caesar.zooniverse.org/workflows/6122/](https://caesar.zooniverse.org/workflows/6122/). 

Replace the workflow number with:
* 6122 for DECALS DR5
* 5653 for SDSS MANGA lost set beta

In the top right, select 'data export' to download.

