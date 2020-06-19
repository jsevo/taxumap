# Common Expected Errors


## Your `taxonomy.csv` file

* Make sure that the `taxonomy.csv` file only contains the columns required for taxumap:

| Kingdom 	| Phylum 	|  ... 	|  Genus 	 	| ASV/OTU 	|
|:-:	|---	|---	|---	|---	
| 'Bacteria" 	| Firmicutes  	| ... 	| Veillonella   	| OTU2352   	| ...   	|

* In other words, make sure your `csv` file does NOT contain an index row

> Often your `taxnomy.csv` dataframe pay contain an index column such as:

|   	| Kingdom  	| Phylum     	| … 	| Genus       	| OTU       	|
|---	|----------	|------------	|---	|-------------	|-----------	|
| 0 	| Bacteria 	| Firmicutes 	| … 	| Veillonella 	| OTU2352 	|
> Notice the  empty column on the left comprising an integer index. This must be removed before calling `taxumap` or `parse_taxonomy_data`. 

> If using Pandas `df.to_csv()` function to create your `csv`, pass it the parameter `index = False` to avoid saving a numerical index.

* Make sure that there are no non-sensical entries, such as:
    1. Mitochondria/Choloplasts
    2. Eukaryotes/Archea
    3. ...