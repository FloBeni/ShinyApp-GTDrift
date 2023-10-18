#You can upload google image of each species with the following bash script


list=$(awk -F '\t' '{print $1}'  ../../../database/list_species.tab)

for query in ${list}; do
echo "${query}"

count=1
find=true

while ${find} = true:
do
find=false



imagelink=$(wget --user-agent 'Mozilla/5.0' -qO - "www.google.be/search?q=${query}\&tbm=isch" | sed 's/</\n</g' | grep '<img' | head -n"$count" | tail -n1 | sed 's/.*src="\([^"]*\)".*/\1/')



if wget $imagelink -O ${query}.png ; then
    echo "Command succeeded"
    wget $imagelink -O ${query}.png
else
    echo "Command failed"
    count=2
    find=true
fi

done
done


