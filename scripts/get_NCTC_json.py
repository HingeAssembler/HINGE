import urllib2
from bs4 import BeautifulSoup
import json

response = urllib2.urlopen('http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/')
html = response.read()

soup=BeautifulSoup(html)
table = soup.find("table")

headings = [th.get_text() for th in table.find("tr").find_all("th")]

dataset={}
for row in table.find_all("tr")[1:]:
    #
    print row 
    row1=  [td.get_text() for td in row.find_all("td")]
    print row1
    metadata={}
    cellname=''
    for i, td in enumerate(row.find_all("td")):
        #print metadata
        link=td.find('a')
        # print i, td
        
        
        if i==1:
            cellname=td.get_text()
            print cellname
        
        if i==3:
            # print td
#             ERR_soup=BeautifulSoup(td)
            ERR_links=[]
            potential_links = td.findAll('a')
            # print potential_links
            for potential_link in td.findAll('a'):
                
                ERR_links.append((potential_link.text, potential_link['href']))
            metadata[headings[i]]=ERR_links
            continue
        if link != None:
            link=link.get('href')
        metadata[headings[i]]=(td.get_text(),link)
    
    list_of_files={}
    for run in metadata[headings[3]]:
        link_to_go=run[1]
        response1 = urllib2.urlopen(link_to_go+"&display=xml")
        xml = response1.read()
        xmlsoup = BeautifulSoup(xml)
        fllist=[]
        for data_block in xmlsoup.findAll('data_block'):
            for files in data_block.findAll('files'):
                for fle in files.findAll('file'):
                    fllist.append(fle['filename'])
        list_of_files[run[0]]=fllist
    metadata['file_paths']=list_of_files
#     print xml
    dataset[cellname]=metadata


with open('NCTC.json', 'w') as outfile:
    json.dump(dataset, outfile)

